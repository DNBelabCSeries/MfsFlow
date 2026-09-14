#!/usr/bin/env python3
"""Small pigz-compatible fallback for hosts with an incompatible libz.

MfsFlow only relies on pigz for stdin/stdout compression and decompression.
This fallback intentionally implements that small interface; it is slower than
pigz and is selected only after the native executable fails a runtime probe.
"""

import gzip
import shutil
import sys
import zlib


def _is_decompression(args):
    return any(
        arg in {"-d", "--decompress"}
        or (arg.startswith("-") and not arg.startswith("--") and "d" in arg[1:])
        for arg in args
    )


def _uses_stdout(args):
    return any(
        arg in {"-c", "--stdout", "--to-stdout"}
        or (arg.startswith("-") and not arg.startswith("--") and "c" in arg[1:])
        for arg in args
    )


def _input_paths(args):
    paths = []
    skip_next = False
    options_ended = False
    for arg in args:
        if skip_next:
            skip_next = False
            continue
        if arg == "--":
            options_ended = True
            continue
        if options_ended:
            paths.append(arg)
            continue
        if arg in {"-p", "--processes", "-b", "--blocksize"}:
            skip_next = True
            continue
        if arg.startswith("--processes=") or arg.startswith("--blocksize="):
            continue
        if arg.startswith("-"):
            continue
        paths.append(arg)
    return paths


def _copy_stream(source, target):
    shutil.copyfileobj(source, target, length=1024 * 1024)


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if any(arg in {"-V", "--version"} for arg in args):
        print("mfsflow pigz compatibility fallback")
        return 0

    decompress = _is_decompression(args)
    to_stdout = _uses_stdout(args) or not _input_paths(args)
    paths = _input_paths(args)

    if decompress:
        if paths:
            for path in paths:
                with gzip.open(path, "rb") as source:
                    _copy_stream(source, sys.stdout.buffer)
        else:
            with gzip.GzipFile(fileobj=sys.stdin.buffer, mode="rb") as source:
                _copy_stream(source, sys.stdout.buffer)
        return 0

    if not to_stdout:
        for path in paths:
            with open(path, "rb") as source, gzip.open(path + ".gz", "wb") as target:
                _copy_stream(source, target)
        return 0

    with gzip.GzipFile(fileobj=sys.stdout.buffer, mode="wb", mtime=0) as target:
        if paths:
            for path in paths:
                with open(path, "rb") as source:
                    _copy_stream(source, target)
        else:
            _copy_stream(sys.stdin.buffer, target)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (BrokenPipeError, OSError, EOFError, gzip.BadGzipFile, zlib.error) as exc:
        print(f"pigz compatibility fallback failed: {exc}", file=sys.stderr)
        raise SystemExit(1)
