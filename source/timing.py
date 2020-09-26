import time
from source import *
from layer import *

testLayer = Layer(4)
Ntests = 100000

testLayer = Layer(1.1)

startTime = time.time()
for i in range(Ntests):
    testLayer = Layer(1 + i)
endTime = time.time()
totalTime = endTime - startTime
timePerTestLayer = totalTime / Ntests

startTime = time.time()
for i in range(Ntests):
    testStack = LayerStack(testLayer, testLayer, testLayer)
endTime = time.time()
totalTime = endTime - startTime
timePerTestStack = totalTime / Ntests

startTime = time.time()
for i in range(Ntests):
    source = Source(wavelength=1, theta=0.5, phi=0.3, pTEM=[1,0.1],layer=testLayer)
endTime = time.time()
totalTime = endTime - startTime
timePerTestSource = totalTime / Ntests

print(f"Total Time per layer creation: {timePerTestLayer}")
print(f"Total Time per stack creation: {timePerTestStack}")
print(f"Total Time per source creation: {timePerTestSource}")
