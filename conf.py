import subprocess, sys, os

DOC_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), 'doc'))
DOC_BUILD = os.path.join(DOC_ROOT, 'doxyxml')

def run_doxygen(folder):
    """Run the doxygen make command in the designated folder"""

    try:
        retcode = subprocess.call("cd %s; make" % folder, shell=True)
        if retcode < 0:
            sys.stderr.write("doxygen terminated by signal %s" % (-retcode))
    except OSError as e:
        sys.stderr.write("doxygen execution failed: %s" % e)


def generate_doxygen_xml(app):
    """Run the doxygen make commands if we're on the ReadTheDocs server"""

    read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

    if read_the_docs_build:

        run_doxygen("../../examples/doxygen")
        run_doxygen("../../examples/specific")
        run_doxygen("../../examples/tinyxml")


def setup(app):

    # Add hook for building doxygen xml when needed
    app.connect("builder-inited", generate_doxygen_xml)
