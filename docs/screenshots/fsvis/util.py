"""
Utility code imported by other modules
"""

def fsvis_main(m, name="screen"):
    from idaes.ui.fsvis import visualize
    import sys, time

    result = visualize(m.fs, name=name, browser=False)
    print(f"[scraper] http://localhost:{result.port}/app?id={name}\n")
    sys.stdout.flush()
    time.sleep(60)
    return 0