g++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) SD3.cpp -o SD3.so
