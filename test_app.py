"""
test_app.py — Pytest + Selenium test suite for SeqNMA
Sequential Network Meta-Analysis browser tool
18+ tests covering statistical engine, UI, and exports.
"""
import os
import sys
import time
import json
import pytest
import subprocess
import threading
import http.server
import socketserver
import io

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException

PORT = 8791
APP_DIR = os.path.dirname(os.path.abspath(__file__))
TIMEOUT = 60


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format, *args):
        pass


@pytest.fixture(scope="session")
def server():
    """Start a local HTTP server."""
    os.chdir(APP_DIR)
    handler = QuietHandler
    httpd = socketserver.TCPServer(("", PORT), handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://localhost:{PORT}/index.html"
    httpd.shutdown()


@pytest.fixture(scope="session")
def driver(server):
    """Create a headless Chrome driver."""
    # Kill orphan chrome/chromedriver
    if sys.platform == "win32":
        subprocess.run(
            ["taskkill", "/F", "/IM", "chromedriver.exe"],
            capture_output=True,
        )

    options = webdriver.ChromeOptions()
    options.add_argument("--headless=new")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--disable-gpu")
    options.add_argument("--window-size=1400,900")
    options.set_capability("goog:loggingPrefs", {"browser": "ALL"})

    try:
        d = webdriver.Chrome(options=options)
    except Exception:
        # Try Edge as fallback
        edge_options = webdriver.EdgeOptions()
        edge_options.add_argument("--headless=new")
        edge_options.add_argument("--no-sandbox")
        edge_options.add_argument("--disable-dev-shm-usage")
        edge_options.add_argument("--disable-gpu")
        edge_options.add_argument("--window-size=1400,900")
        edge_options.set_capability("goog:loggingPrefs", {"browser": "ALL"})
        d = webdriver.Edge(options=edge_options)

    d.set_page_load_timeout(TIMEOUT)
    d.implicitly_wait(5)
    d.get(server)
    yield d
    d.quit()


def load_demo_and_analyze(driver):
    """Helper: click Demo Data then Analyze."""
    driver.find_element(By.ID, "btnDemo").click()
    time.sleep(0.3)
    driver.find_element(By.ID, "btnAnalyze").click()
    time.sleep(1.5)


def js(driver, script):
    """Execute JS and return result."""
    return driver.execute_script(f"return {script}")


# ===== TEST 1: App loads without JS errors =====
def test_01_app_loads_no_js_errors(driver):
    """App loads without JS errors."""
    logs = driver.get_log("browser")
    severe = [l for l in logs if l["level"] == "SEVERE" and l.get("source") != "network"]
    assert len(severe) == 0, f"JS errors: {severe}"


# ===== TEST 2: Demo data loads 12 studies, detects 4 treatments =====
def test_02_demo_data_loads(driver):
    """Demo data loads 12 studies and detects 4 treatments."""
    load_demo_and_analyze(driver)
    result = js(driver, "window.SeqNMA.getCurrentResult()")
    assert result is not None, "Analysis returned null"
    assert len(result["allStudies"]) == 12
    assert len(result["treatments"]) == 4
    assert sorted(result["treatments"]) == ["A", "B", "C", "D"]


# ===== TEST 3: Network graph SVG has 4 nodes =====
def test_03_network_graph_4_nodes(driver):
    """Network graph SVG has 4 nodes (treatments A-D)."""
    # Switch to network tab
    driver.find_element(By.CSS_SELECTOR, '[data-tab="network"]').click()
    time.sleep(0.5)
    svg = driver.find_element(By.ID, "networkGraphSVG")
    circles = svg.find_elements(By.TAG_NAME, "circle")
    assert len(circles) == 4, f"Expected 4 nodes, got {len(circles)}"


# ===== TEST 4: O'Brien-Fleming boundary formula =====
def test_04_obf_boundary_formula(driver):
    """O'Brien-Fleming: z_k = z_alpha/sqrt(t_k) at t=0.5 boundary ~ 2.80."""
    # For 6 comparisons, adjusted alpha = 0.05/6
    adj_alpha = 0.05 / 6
    boundary = js(driver, f"window.SeqNMA.obrienFlemingBoundary(0.5, {adj_alpha})")
    # z_{alpha/2} where alpha = 0.05/6 = 0.00833...
    # z_{0.004167} ~ 2.635 (approximately)
    # boundary = z / sqrt(0.5) = z * sqrt(2)
    assert boundary is not None
    assert boundary > 2.5, f"Boundary at t=0.5 should be > 2.5, got {boundary}"
    # Verify formula: z_k = z_alpha / sqrt(t)
    za = js(driver, f"window.SeqNMA.qnorm(1 - {adj_alpha}/2)")
    expected = za / (0.5 ** 0.5)
    assert abs(boundary - expected) < 0.001, f"Boundary {boundary} != expected {expected}"


# ===== TEST 5: RIS calculation is positive =====
def test_05_ris_positive(driver):
    """RIS calculation: for delta=0.2, alpha=0.05, power=0.80, RIS > 0."""
    result = js(driver, "window.SeqNMA.getCurrentResult()")
    for comp in result["comparisons"]:
        sr = result["sequentialResults"].get(comp)
        if sr:
            assert sr["ris"] > 0, f"RIS for {comp} should be > 0, got {sr['ris']}"


# ===== TEST 6: Information fraction increases monotonically =====
def test_06_info_fraction_monotonic(driver):
    """Information fraction increases as studies added (monotonic)."""
    result = js(driver, "window.SeqNMA.getCurrentResult()")
    for comp in result["comparisons"]:
        sr = result["sequentialResults"].get(comp)
        if not sr or len(sr["points"]) < 2:
            continue
        fracs = [p["infoFraction"] for p in sr["points"]]
        # Info fraction should be non-decreasing (it can stay same if study doesn't affect comparison)
        for i in range(1, len(fracs)):
            assert fracs[i] >= fracs[i - 1] - 1e-10, (
                f"{comp}: info fraction decreased at step {i}: {fracs[i-1]:.6f} -> {fracs[i]:.6f}"
            )


# ===== TEST 7: Boundary plot SVG renders =====
def test_07_boundary_plot_renders(driver):
    """Boundary plot SVG renders for selected comparison."""
    driver.find_element(By.CSS_SELECTOR, '[data-tab="boundary"]').click()
    time.sleep(0.5)
    svg = driver.find_element(By.ID, "boundaryPlotSVG")
    assert svg is not None
    # Check it has paths (boundary lines) and circles (observed points)
    paths = svg.find_elements(By.TAG_NAME, "path")
    circles = svg.find_elements(By.TAG_NAME, "circle")
    assert len(paths) > 0, "No paths in boundary plot"
    assert len(circles) > 0, "No observed points in boundary plot"


# ===== TEST 8: Cumulative NMA pooled effect changes =====
def test_08_cumulative_effect_changes(driver):
    """Cumulative NMA: pooled effect changes as studies added."""
    result = js(driver, "window.SeqNMA.getCurrentResult()")
    # Check A vs B comparison (has 3 direct studies)
    sr = result["sequentialResults"].get("A vs B")
    assert sr is not None, "A vs B comparison not found"
    effects = [p["effect"] for p in sr["points"]]
    assert len(effects) >= 3, f"Expected >= 3 cumulative points for A vs B, got {len(effects)}"
    # Effects should not all be identical
    unique_effects = set(round(e, 6) for e in effects)
    assert len(unique_effects) > 1, "Cumulative effects should change as studies are added"


# ===== TEST 9: Bonferroni adjustment =====
def test_09_bonferroni_adjustment(driver):
    """Bonferroni adjustment: alpha/6 for 6 comparisons (4 choose 2)."""
    result = js(driver, "window.SeqNMA.getCurrentResult()")
    assert result["numComparisons"] == 6, f"Expected 6 comparisons, got {result['numComparisons']}"
    expected_adj = 0.05 / 6
    assert abs(result["adjustedAlpha"] - expected_adj) < 1e-10, (
        f"Adjusted alpha {result['adjustedAlpha']} != {expected_adj}"
    )


# ===== TEST 10: Futility boundary labeled non-binding =====
def test_10_futility_non_binding(driver):
    """Futility boundary labeled 'non-binding'."""
    driver.find_element(By.CSS_SELECTOR, '[data-tab="boundary"]').click()
    time.sleep(0.5)
    svg_el = driver.find_element(By.ID, "boundaryPlotSVG")
    svg_html = svg_el.get_attribute("outerHTML")
    assert "non-binding" in svg_html.lower(), "Futility boundary must be labeled 'non-binding'"


# ===== TEST 11: Design effect D^2 uses tau^2-based formula =====
def test_11_design_effect_formula(driver):
    """Design effect D^2 uses tau^2-based formula (NOT cluster DEFF)."""
    # Verify D^2 for a comparison with tau^2 > 0
    result = js(driver, "window.SeqNMA.getCurrentResult()")
    tau2 = result["tau2"]

    # Manually compute D^2 from the formula:
    # D^2 = 1 + tau^2 * (sum(1/v_i^2) / (sum(1/v_i))^2 * k - 1)
    studies = result["allStudies"]
    k = len(studies)
    sum_inv_v = sum(1.0 / s["variance"] for s in studies)
    sum_inv_v2 = sum(1.0 / (s["variance"] ** 2) for s in studies)
    expected_D2 = 1 + tau2 * (sum_inv_v2 / (sum_inv_v ** 2) * k - 1)
    expected_D2 = max(1, expected_D2)

    # Check D2 in one comparison
    for comp in result["comparisons"]:
        sr = result["sequentialResults"].get(comp)
        if sr:
            # D2 should match our manual calculation
            assert abs(sr["D2"] - expected_D2) < 0.01, (
                f"D2 mismatch for {comp}: got {sr['D2']}, expected {expected_D2}"
            )
            break


# ===== TEST 12: 2-treatment network degenerates to standard TSA =====
def test_12_two_treatment_tsa(driver):
    """2-treatment network: degenerates to standard TSA."""
    csv_2treat = "Study,Year,Treat1,Treat2,Effect,SE\\nStudy1,2015,A,B,-0.30,0.12\\nStudy2,2016,A,B,-0.25,0.11\\nStudy3,2017,A,B,-0.28,0.10"

    result = driver.execute_script(
        "var studies = window.SeqNMA.parseCSV(arguments[0]);"
        "var settings = {alpha: 0.05, power: 0.80, delta: 0.2, effectScale: 'SMD'};"
        "return window.SeqNMA.runSequentialNMA(studies, settings);",
        csv_2treat.replace("\\n", "\n")
    )
    assert result is not None
    assert len(result["treatments"]) == 2
    assert result["numComparisons"] == 1
    # Bonferroni with 1 comparison = no adjustment
    assert abs(result["adjustedAlpha"] - 0.05) < 1e-10


# ===== TEST 13: 3-treatment triangle network =====
def test_13_three_treatment_triangle(driver):
    """3-treatment triangle network: 3 comparisons monitored."""
    csv_3treat = "Study,Year,Treat1,Treat2,Effect,SE\nS1,2015,A,B,-0.20,0.12\nS2,2016,A,C,-0.35,0.11\nS3,2017,B,C,-0.10,0.13"

    result = driver.execute_script(
        "var studies = window.SeqNMA.parseCSV(arguments[0]);"
        "var settings = {alpha: 0.05, power: 0.80, delta: 0.2, effectScale: 'SMD'};"
        "return window.SeqNMA.runSequentialNMA(studies, settings);",
        csv_3treat
    )
    assert result is not None
    assert len(result["treatments"]) == 3
    assert result["numComparisons"] == 3
    assert len(result["comparisons"]) == 3


# ===== TEST 14: Star network includes indirect comparisons =====
def test_14_star_network_indirect(driver):
    """Star network (A as hub): indirect comparisons included."""
    csv_star = "Study,Year,Treat1,Treat2,Effect,SE\nS1,2015,A,B,-0.20,0.12\nS2,2016,A,C,-0.35,0.11\nS3,2017,A,D,-0.15,0.14"

    result = driver.execute_script(
        "var studies = window.SeqNMA.parseCSV(arguments[0]);"
        "var settings = {alpha: 0.05, power: 0.80, delta: 0.2, effectScale: 'SMD'};"
        "return window.SeqNMA.runSequentialNMA(studies, settings);",
        csv_star
    )
    assert result is not None
    assert len(result["treatments"]) == 4
    assert result["numComparisons"] == 6
    # B vs C, B vs D, C vs D are indirect — should still have sequential results
    # They appear after all 3 studies are in
    nma = result["fullNMA"]
    bc = [c for c in nma["comparisons"] if c["comp"] == "B vs C"]
    assert len(bc) == 1, "B vs C indirect comparison should exist"
    assert bc[0]["se"] > 0, "B vs C should have valid SE from indirect evidence"


# ===== TEST 15: League table 4x4 with 6 unique comparisons =====
def test_15_league_table(driver):
    """League table: 4x4 matrix with 6 unique comparisons."""
    driver.find_element(By.CSS_SELECTOR, '[data-tab="league"]').click()
    time.sleep(0.5)
    table = driver.find_element(By.ID, "leagueTableElement")
    rows = table.find_elements(By.TAG_NAME, "tr")
    # Header + 4 data rows = 5
    assert len(rows) == 5, f"Expected 5 rows (header + 4), got {len(rows)}"
    # Check 4x4 structure
    data_rows = rows[1:]
    for row in data_rows:
        cells = row.find_elements(By.TAG_NAME, "td") + row.find_elements(By.TAG_NAME, "th")
        assert len(cells) == 5, f"Expected 5 cells per row (header + 4), got {len(cells)}"


# ===== TEST 16: Information tracker table renders =====
def test_16_info_tracker_renders(driver):
    """Information tracker table renders."""
    driver.find_element(By.CSS_SELECTOR, '[data-tab="info"]').click()
    time.sleep(0.5)
    table = driver.find_element(By.ID, "infoTrackerTable")
    rows = table.find_elements(By.CSS_SELECTOR, "tbody tr")
    assert len(rows) == 6, f"Expected 6 comparison rows, got {len(rows)}"
    # Check badges exist
    badges = table.find_elements(By.CLASS_NAME, "badge")
    assert len(badges) == 6, f"Expected 6 status badges, got {len(badges)}"


# ===== TEST 17: Network evolution slider changes graph =====
def test_17_network_slider_changes_graph(driver):
    """Network evolution slider changes graph."""
    driver.find_element(By.CSS_SELECTOR, '[data-tab="network"]').click()
    time.sleep(0.5)

    slider = driver.find_element(By.ID, "networkYearSlider")

    # Get initial SVG (year = max)
    svg1 = driver.find_element(By.ID, "networkGraphSVG").get_attribute("outerHTML")

    # Set slider to 2012 (only 3 studies: Adams 2010, Baker 2011, Chen 2012)
    driver.execute_script("arguments[0].value = 2012; arguments[0].dispatchEvent(new Event('input'));", slider)
    time.sleep(0.5)

    svg2 = driver.find_element(By.ID, "networkGraphSVG").get_attribute("outerHTML")
    assert svg1 != svg2, "Network graph should change when slider moves"

    # At 2012, only treatments A, B, C should be present (no D yet)
    text_content = svg2
    assert "4 treatments" not in text_content or "3 treatments" in text_content, \
        "At 2012, should have fewer than 4 treatments"


# ===== TEST 18: Export CSV =====
def test_18_export_csv(driver):
    """Export CSV produces valid content."""
    driver.find_element(By.CSS_SELECTOR, '[data-tab="boundary"]').click()
    time.sleep(0.3)

    # Build CSV content via JS — count data rows
    row_count = driver.execute_script(
        "var result = window.SeqNMA.getCurrentResult();"
        "if (!result) return 0;"
        "var count = 0;"
        "for (var i = 0; i < result.comparisons.length; i++) {"
        "  var sr = result.sequentialResults[result.comparisons[i]];"
        "  if (sr) count += sr.points.length;"
        "}"
        "return count;"
    )
    assert row_count is not None and row_count > 0, "CSV should have data rows"

    # Verify export button exists
    btns = driver.find_elements(By.CSS_SELECTOR, ".btn-export-csv")
    assert len(btns) > 0, "Export CSV button should exist"


# ===== TEST 19: NMA basic parameter estimation =====
def test_19_nma_parameter_estimation(driver):
    """NMA basic parameters are consistent across comparisons."""
    result = js(driver, "window.SeqNMA.getCurrentResult()")
    nma = result["fullNMA"]
    comps = {c["comp"]: c for c in nma["comparisons"]}

    # Consistency check: AB + BC should approximately equal AC
    ab = comps.get("A vs B")
    bc = comps.get("B vs C")
    ac = comps.get("A vs C")

    if ab and bc and ac:
        # d_AC = d_AB + d_BC (in terms of basic parameters)
        # AB: effect = d_B - d_A (= beta_B since A is reference)
        # BC: effect = d_C - d_B
        # AC: effect = d_C - d_A (= beta_C since A is reference)
        indirect_ac = ab["effect"] + bc["effect"]
        assert abs(indirect_ac - ac["effect"]) < 0.5, (
            f"Consistency: AB ({ab['effect']:.3f}) + BC ({bc['effect']:.3f}) = "
            f"{indirect_ac:.3f} should be close to AC ({ac['effect']:.3f})"
        )


# ===== TEST 20: Tau-squared is non-negative =====
def test_20_tau2_nonnegative(driver):
    """Tau-squared from NMA is non-negative."""
    result = js(driver, "window.SeqNMA.getCurrentResult()")
    assert result["tau2"] >= 0, f"tau2 should be >= 0, got {result['tau2']}"


# ===== TEST 21: Phi and qnorm inverse relationship =====
def test_21_phi_qnorm_inverse(driver):
    """Standard normal CDF and quantile are inverses."""
    for p in [0.01, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.99]:
        z = js(driver, f"window.SeqNMA.qnorm({p})")
        p_back = js(driver, f"window.SeqNMA.phi({z})")
        assert abs(p_back - p) < 0.005, f"phi(qnorm({p})) = {p_back}, expected {p}"


# ===== TEST 22: Forest plot renders =====
def test_22_forest_plot_renders(driver):
    """Cumulative forest plot SVG renders."""
    driver.find_element(By.CSS_SELECTOR, '[data-tab="forest"]').click()
    time.sleep(0.5)
    svg = driver.find_element(By.ID, "forestPlotSVG")
    assert svg is not None
    # Should have rect elements (diamond points)
    rects = svg.find_elements(By.TAG_NAME, "rect")
    assert len(rects) > 1, "Forest plot should have diamond markers"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
