if ! test -d build; then mkdir build; fi
cd build
cmake .. && make
cd ..
python setup.py install
