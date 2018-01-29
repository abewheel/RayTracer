# Ray Tracer


<p align="center">

  <img src="https://github.com/abewheel/RayTracer/assign3/blob/master/003.jpg" alt="Shelves"/>

</p>



## Description


Ray Tracer is a Visual Studio 2017 C++ program that accepts a .txt file defining a 3D object space and outputs an image of the space.


<p align="center">

  <img src="https://github.com/abewheel/RayTracer/assign3/blob/master/004.jpg" alt="Text"/>

</p>



## About the Project


This repository contains the Visual Studio 2017 project code for Ray Tracer. I wrote Ray Tracer in 2017 at USC as a part of my Computer Graphics course with Dr. Hao Li, CSCI420.


<p align="center">

  <img src="https://github.com/abewheel/RayTracer/assign3/blob/master/001.jpg" alt="Ball and triangle"/>

</p>



## Technical Details


Roller Coaster was implemented in Visual Studio 2017 using C++ and OpenGL. Most of the action occurs in RayTracer/assign3/assign3.cpp. It accepts a list of objects defined as either spheres, triangles, or lights. Each object is accompanied with a description of its position, diffuse & specular color, shininess coefficient, radius, and normal (each where applicable) and utilizes the Phong shading technique to determine the RGB values for each pixel. 