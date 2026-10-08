"""The PhLynx manifests: now cam_testing/manifests.py (python -m cam_testing.manifests); kept for old commands."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cam_testing.manifests import main  # noqa: E402

if __name__ == '__main__':
    main()
