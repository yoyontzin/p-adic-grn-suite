from padic_grn_viz import __version__


def test_version_defined() -> None:
    assert isinstance(__version__, str)
    assert __version__ == "0.1.0"
