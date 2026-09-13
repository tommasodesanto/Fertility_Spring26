# September slides: automatic PDF rebuild

`watch_september_slides.py` runs latexmk's continuous dependency watcher on
`latex/september_14_presentation.tex`. Saving the source or a detected input in
any editor rebuilds the deck. Successful builds atomically replace each PDF at
`latex/september_14_presentation.pdf` and `output/pdf/september_14_presentation.pdf`.
Failed builds keep the last successfully published PDFs. Viewer refresh depends
on the PDF viewer; reopen the same PDF if it does not reload changed files.

Auxiliary files and logs stay in `tmp/september_slides_autobuild/`. Latexmk makes
as many passes as references require. Its local `-norc` option avoids the
repository's older in-place output setting. No other deck is watched.

The watcher runs as a detached process started from the application session,
with its PID in `tmp/september_slides_autobuild/watcher.pid` and its log in
`watcher.log` alongside it. macOS denied a login LaunchAgent access to the
Desktop project; that unsuccessful service registration was removed. This is
not configured to restart after a logout or reboot. There should be only one
watcher. While it is active,
editors should let it publish the shared PDFs rather than copy an older build
over them. Manual verification can use a separate scratch output directory.

Stop:

```sh
kill "$(cat tmp/september_slides_autobuild/watcher.pid)"
```

Start again:

```sh
mkdir -p tmp/september_slides_autobuild
nohup python3 ops/latex-preview/watch_september_slides.py > tmp/september_slides_autobuild/watcher.log 2>&1 < /dev/null &
echo $! > tmp/september_slides_autobuild/watcher.pid
```

Run these commands from the project root. The tracked Python script can also
be run in the foreground from a terminal. Check for an existing watcher before
starting another one.
