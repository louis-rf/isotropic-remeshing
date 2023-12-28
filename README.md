
I installed geogram from here: https://github.com/BrunoLevy/geogram/wiki

For OSX this is:
```
git clone --recurse-submodules https://github.com/BrunoLevy/geogram.git
cd geogram
./configure.sh
cd Build/Darwin-clang-dynamic-Release
make -j 8
```

I based the pygeogram package off the method shown [here](https://github.com/pybind/cmake_example) to make python bindings for
c++.
```
git clone --recursive https://github.com/pybind/cmake_example.git
# move the required files, i.e. all the files which aren't already there.
for item in cmake_example/*; do [ ! -e "pygeogram/$(basename "$item")" ] && mv "$item" "pygeogram/"; done
rm -rf cmake_example
```

Now build the binding:
```
cd pygeogram
pip install .
```