from pathlib import Path

from playwright.sync_api import sync_playwright


def main() -> None:
    root = Path(__file__).resolve().parent
    frontend_file = root / "frontend" / "index.html"
    output_dir = root / "screenshots"
    output_dir.mkdir(parents=True, exist_ok=True)

    url = frontend_file.as_uri()

    with sync_playwright() as p:
        browser = p.chromium.launch()

        page = browser.new_page(viewport={"width": 1440, "height": 900})
        page.goto(url, wait_until="networkidle")
        page.screenshot(path=str(output_dir / "demo-ui-desktop.png"), full_page=True)
        page.close()

        page_mobile = browser.new_page(
            viewport={"width": 430, "height": 932},
            device_scale_factor=2,
            is_mobile=True,
            has_touch=True,
        )
        page_mobile.goto(url, wait_until="networkidle")
        page_mobile.screenshot(path=str(output_dir / "demo-ui-mobile.png"), full_page=True)
        page_mobile.close()

        browser.close()


if __name__ == "__main__":
    main()
