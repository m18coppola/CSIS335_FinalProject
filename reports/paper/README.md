# Ray Tracing Parallelization

Use this file to keep track of your thoughts and efforts!
The more we can document here, the less we have to use our brains when we have to regurgitate it into a LaTex document.

## Interesting Links
- tinyraytracer
	- https://github.com/ssloy/tinyraytracer/wiki/
- Portable Pixmap Format
	- https://en.wikipedia.org/wiki/Netpbm
- cglm (header-only graphics math lib for c)
	- https://cglm.readthedocs.io/en/latest/
- Random distribution of rays among threads
	- https://ieeexplore.ieee.org/abstract/document/1580017
- General raytracing optimizations
	- https://ieeexplore.ieee.org/abstract/document/537389
- Distribution of data for complex scenes
	- https://ieeexplore.ieee.org/abstract/document/291533
- KD-tree data divisions
	- https://arthurflor23.medium.com/ray-tracing-with-kd-tree-from-scratch-9b218823bc00
- triangle-ray intersection
	- https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution

## Further Research:
- extreme parallelization by minimizing mutex locks with complex jobs (bag of bags of tasks?)
  - Faster with predefining task trees to reduce mutex locks?
  - Split tree when worker is caught idling?
- processing chunks of image per thread (seems to be used in industry for some reason)


