"Various utilities used by sourmash tests."

import sys
import os
import tempfile
import shutil
import subprocess
import collections
import pprint

import importlib.metadata
import traceback
from io import open  # pylint: disable=redefined-builtin
from io import StringIO


def get_test_data(filename):
    thisdir = os.path.dirname(__file__)
    return os.path.join(thisdir, "test-data", filename)


def make_file_list(filename, paths):
    with open(filename, "wt") as fp:
        fp.write("\n".join(paths))
        fp.write("\n")


def zip_siglist(runtmp, siglist, db):
    runtmp.sourmash("sig", "cat", siglist, "-o", db)
    return db


def index_siglist(
    runtmp,
    siglist,
    db,
    *,
    ksize=31,
    scaled=None,
    moltype="DNA",
    toggle_internal_storage="--internal-storage",
):
    # build index
    extra_args = []
    if scaled is not None:
        extra_args = ["--scaled", str(scaled)]

    runtmp.sourmash(
        "scripts",
        "index",
        siglist,
        "-o",
        db,
        "-k",
        str(ksize),
        "--moltype",
        moltype,
        toggle_internal_storage,
        *extra_args,
    )
    return db


def _runscript(scriptname):
    """Find & run a script with exec (i.e. not via os.system or subprocess)."""
    namespace = {"__name__": "__main__"}
    namespace["sys"] = globals()["sys"]

    entry_points = importlib.metadata.entry_points(
        group="console_scripts", name="sourmash"
    )
    assert len(entry_points) == 1
    smash_cli = tuple(entry_points)[0].load()
    smash_cli()
    return 0


ScriptResults = collections.namedtuple("ScriptResults", ["status", "out", "err"])


def runscript(scriptname, args, **kwargs):
    """Run a Python script using exec().

    Run the given Python script, with the given args, in the given directory,
    using 'exec'.  Mimic proper shell functionality with argv, and capture
    stdout and stderr.

    When using :attr:`fail_ok`=False in tests, specify the expected error.
    """
    __tracebackhide__ = True
    sysargs = [scriptname]
    sysargs.extend(args)

    cwd = os.getcwd()
    in_directory = kwargs.get("in_directory", cwd)
    fail_ok = kwargs.get("fail_ok", False)

    try:
        status = -1
        oldargs = sys.argv
        sys.argv = sysargs

        oldin = None
        if "stdin_data" in kwargs:
            oldin, sys.stdin = sys.stdin, StringIO(kwargs["stdin_data"])

        oldout, olderr = sys.stdout, sys.stderr
        sys.stdout = StringIO()
        sys.stdout.name = "StringIO"
        sys.stderr = StringIO()

        os.chdir(in_directory)

        try:
            print("running:", scriptname, "in:", in_directory, file=oldout)
            print("arguments", sysargs, file=oldout)

            status = _runscript(scriptname)
        except SystemExit as err:
            status = err.code
            if status == None:
                status = 0
        except:  # pylint: disable=bare-except
            traceback.print_exc(file=sys.stderr)
            status = -1
    finally:
        sys.argv = oldargs
        out, err = sys.stdout.getvalue(), sys.stderr.getvalue()
        sys.stdout, sys.stderr = oldout, olderr

        if oldin:
            sys.stdin = oldin

        os.chdir(cwd)

    if status != 0 and not fail_ok:
        print(out)
        print(err)
        assert False, (status, out, err)

    return ScriptResults(status, out, err)


class TempDirectory(object):
    def __init__(self):
        self.tempdir = tempfile.mkdtemp(prefix="sourmashtest_")

    def __enter__(self):
        return self.tempdir

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            shutil.rmtree(self.tempdir, ignore_errors=True)
        except OSError:
            pass

        if exc_type:
            return False


class SourmashCommandFailed(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)
        self.message = msg


class RunnerContext(object):
    """
    I am a RunnerContext object from sourmash_tst_utils.

    I have methods 'run_sourmash' and 'run', which run Python scripts.

    Take a look at my 'location', 'last_command' and 'last_result' attributes!

    You can use the 'output' method to build filenames in my temp directory.
    """

    def __init__(self, location):
        self.location = location
        self.last_command = None
        self.last_result = None

    def run_sourmash(self, *args, **kwargs):
        "Run the sourmash script with the given arguments."
        kwargs["fail_ok"] = True
        if "in_directory" not in kwargs:
            kwargs["in_directory"] = self.location

        cmdlist = ["sourmash"]
        cmdlist.extend((str(x) for x in args))
        self.last_command = " ".join(cmdlist)
        self.last_result = runscript("sourmash", args, **kwargs)

        if self.last_result.status:
            raise SourmashCommandFailed(self.last_result.err)

        return self.last_result

    sourmash = run_sourmash

    def run(self, scriptname, *args, **kwargs):
        "Run a script with the given arguments."
        if "in_directory" not in kwargs:
            kwargs["in_directory"] = self.location
        self.last_command = " ".join(args)
        self.last_result = runscript(scriptname, args, **kwargs)
        return self.last_result

    def output(self, path):
        return os.path.join(self.location, path)

    def __str__(self):
        s = ""
        if self.last_command:
            s += "Last command run:\n{}\n".format(repr(self.last_command))
            if self.last_result:
                s += "\nLAST RESULT:\n"
                s += "- exit code: {}\n\n".format(self.last_result.status)
                if self.last_result.out:
                    s += "- stdout:\n---\n{}---\n".format(self.last_result.out)
                else:
                    s += "(no stdout)\n\n"
                if self.last_result.err:
                    s += "- stderr:\n---\n{}---\n".format(self.last_result.err)
                else:
                    s += "(no stderr)\n"

        return s


def in_tempdir(fn):
    def wrapper(*args, **kwargs):
        with TempDirectory() as location:
            ctxt = RunnerContext(location)
            newargs = [ctxt] + list(args)
            return fn(*newargs, **kwargs)

    return wrapper


def in_thisdir(fn):
    def wrapper(*args, **kwargs):
        ctxt = RunnerContext(os.getcwd())
        newargs = [ctxt] + list(args)
        return fn(*newargs, **kwargs)

    return wrapper
