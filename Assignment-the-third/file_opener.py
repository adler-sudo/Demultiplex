#!/usr/bin/env python

from contextlib import ExitStack


filenames = ["test.txt", "another.txt", "donkey.txt"]

for f in filenames:
    with open(f, 'w') as f:
        f.close()

with ExitStack() as stack:
    files = [stack.enter_context(open(fname)) for fname in filenames]