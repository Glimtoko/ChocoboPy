SO_NAME = lagrangian_hydro.cpython-37m-x86_64-linux-gnu.so
INSTALL_DIR = /prod/chocobopy/bin/
BIN_DIR = /prod/bin


.PHONY: exe install clean fullclean

PYFILES = chocobopy.py \
	geometry.py \
	data_structures.py \
	initial_conditions.py \
	output.py \
	eos.py \
	inputs.py \
	read_mesh.py

$(SO_NAME): lagrangian_hydro.pyx
	python setup.py build_ext --inplace

exe: dist/chocobopy

dist/chocobopy: $(SO_NAME) $(PYFILES)
	pyinstaller --onefile chocobopy.py

install: dist/chocobopy
	cp dist/chocobopy $(INSTALL_DIR)
	chmod +x $(INSTALL_DIR)/chocobopy
	ln -sf $(INSTALL_DIR)/chocobopy $(BIN_DIR)/chocobopy

clean:
	-rm chocobopy.spec
	-rm -r build
	-rm -r __pycache__

fullclean:
	-rm chocobopy.spec
	-rm -r build
	-rm -r __pycache__
	-rm -r dist
	-rm $(SO_NAME)
