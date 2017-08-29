Performance Programming
===

The aim of this coursework is to perform single thread performance optimisation on the back end compute nodes of ARCHER. 

####Program details:

If two particles approach closer than `Size` we flip the direction of the
interaction force to approximate a collision. Coordinates are relative to a large central mass and the entire system is moving relative to the
viscous media. If two particles approach closer than `Size` we flip the direction of the
interaction force to approximate a collision.

The program reads `input.dat` and writes `output.dat`.

The directory Test contains source of a program to compare `output.dat` files.
It will only report differences above a preset tolerance value.





