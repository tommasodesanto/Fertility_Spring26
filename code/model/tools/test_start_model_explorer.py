"""Cheap launcher identity checks; no solution loading or model imports."""
import hashlib
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import socket
import subprocess
import sys
import tempfile
import threading
import unittest


LAUNCHER = Path(__file__).with_name('start_model_explorer.command')
PROBE = LAUNCHER.read_text().split("<<'PY'\n", 1)[1].split('\nPY\n', 1)[0]


class LauncherIdentityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.config = self.root / 'cases.json'
        self.config.write_text('{"cases": [{"id": "latest"}]}')
        # Let the OS select one test port; reserve its next two while checking.
        for _ in range(100):
            sockets = []
            try:
                first = socket.socket()
                sockets.append(first)
                first.bind(('127.0.0.1', 0))
                self.port = first.getsockname()[1]
                for port in (self.port + 1, self.port + 2):
                    sock = socket.socket()
                    sockets.append(sock)
                    sock.bind(('127.0.0.1', port))
                break
            except OSError:
                continue
            finally:
                for sock in sockets:
                    sock.close()
        else:
            self.fail('Could not find three local test ports')

    def serve(self, port, meta):
        class Handler(BaseHTTPRequestHandler):
            def do_GET(self):
                body = json.dumps(meta).encode()
                self.send_response(200)
                self.end_headers()
                self.wfile.write(body)

            def log_message(self, *_):
                pass

        server = ThreadingHTTPServer(('127.0.0.1', port), Handler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()

        def stop():
            server.shutdown()
            server.server_close()
            thread.join()

        self.addCleanup(stop)

    def identity(self, config=None):
        config = (config or self.config).resolve()
        return dict(config_path=str(config),
                    config_sha256=hashlib.sha256(config.read_bytes()).hexdigest())

    def select(self, config=None):
        probe = PROBE.replace('range(8765, 8786)',
                              f'range({self.port}, {self.port + 3})')
        result = subprocess.run([sys.executable, '-', str(config or self.config)],
                                input=probe, text=True, capture_output=True,
                                check=True, timeout=5)
        return result.stdout.splitlines()

    def test_reuses_only_matching_identity(self):
        self.serve(self.port, self.identity())
        self.assertEqual(self.select(), [str(self.config), str(self.port), 'reuse'])

    def test_historical_server_without_identity_is_preserved(self):
        self.serve(self.port, {'cases': [{'id': 'soft'}]})
        expected = [str(self.config), str(self.port + 1), 'start']
        self.assertEqual(self.select(), expected)
        # A second probe proves the existing process was left serving.
        self.assertEqual(self.select(), expected)

    def test_reuses_matching_server_after_vacant_port(self):
        self.serve(self.port + 1, self.identity())
        self.assertEqual(self.select(), [str(self.config), str(self.port + 1), 'reuse'])

    def test_same_path_changed_config_is_not_reused(self):
        self.serve(self.port, self.identity())
        self.config.write_text('{"cases": [{"id": "new"}]}')
        self.assertEqual(self.select(), [str(self.config), str(self.port + 1), 'start'])

    def test_latest_symlink_change_is_not_reused(self):
        latest = self.root / 'latest.json'
        latest.symlink_to(self.config)
        self.serve(self.port, self.identity(latest))
        self.assertEqual(self.select(latest)[2], 'reuse')
        other = self.root / 'new_cases.json'
        other.write_text(self.config.read_text())
        latest.unlink()
        latest.symlink_to(other)
        self.assertEqual(self.select(latest), [str(other), str(self.port + 1), 'start'])


if __name__ == '__main__':
    unittest.main()
