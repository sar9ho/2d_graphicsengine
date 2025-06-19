# 2d_graphicsengine
custom 2D graphics engine with a focus on low-level rendering, color blending, and geometry processing.

To run/test:
    > make
    > ./image -e expected
    > ./tests
    > ./bench
    > ./dbench
    
For more details when testing...
    > mkdir diff
    > ./image -e expected -d diff
    > open diff/index.html
    >
    > ./tests -v
    >
    > ./bench -m substring-to-match
