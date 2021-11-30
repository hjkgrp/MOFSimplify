*Use this info once inside the virtual machine*

To start MOFSimplify, first enter the following command into the terminal:
conda activate MOFSimplify

Then enter the mofSimplify folder in the terminal. Enter the following commands into the terminal:
export PYTHONPATH=${PYTHONPATH}:${PWD}
twistd -n web --port tcp:8000 --wsgi app.app

Note: may need to kill a previous twistd process by using the PID