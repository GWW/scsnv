#!/usr/bin/env python
#Copyright (c) 2018-2020 Gavin W. Wilson
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

def go():
    import sys
    from scsnvpy import cmds

    if len(sys.argv) < 2:
        cmd = ''
    else:
        cmd = sys.argv[1]
        sys.argv[1] = f'{sys.argv[0]} {sys.argv[1]}'
        args = list(sys.argv[1:])

    if cmd not in cmds.COMMANDS:
        print("Invalid program available options are {}".format(", ".join(cmds.COMMANDS.keys())))
        sys.exit()

    cmds.COMMANDS[cmd](args)
