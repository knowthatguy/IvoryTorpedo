import os
import sys
import unittest

import app

cmd_folder = os.path.abspath(os.path.join(".."))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)


class TestIvoryTorpedo(unittest.TestCase):
    def setUp(self):
        app.app.config["TESTING"] = True
        self.app = app.app.test_client()

    def test_main_form(self):
        rv = self.app.get("/")
        self.assertTrue(len(rv.data) > 0)

    def callculate(self):
        return self.app.post("/charts", data=dict(R=0.5, H=1, follow_redirects=True))

    def test_charts(self):
        rv = self.callculate()
        self.assertTrue(len(rv.data) > 0)


if __name__ == "__main__":
    unittest.main()
