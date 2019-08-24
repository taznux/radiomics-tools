if ! test -d build; then mkdir build; fi
cd build
cmake .. && make
cd ..
python3 setup.py install
