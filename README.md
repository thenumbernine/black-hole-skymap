# A NOTE ABOUT ANYONE FOLLOWING THE LINK FROM http://spiro.fisica.unipd.it/~antonell/schwarzschild/

This guy was kind enough to reference this page and publicise a previous math error I had made (I had forgot the cross terms in the Cartesian form of the Schwarzschild metric). If only he was kind enough to attempt to contact me in private first, I could have reminded him that the renderer implementation's geodesic calculations were based on the correct pseudo-Cartesian metric all along. Here's where I got the metric from:  [Catalogue of Spacetimes](https://arxiv.org/abs/0904.4184) by Mueller & Grave. Here's my [geodesic calculations](http://christopheremoore.net/symbolic-lua/test-output/Schwarzschild%20-%20Cartesian.html) based on my [Computer Algebra System](http://christopheremoore.net/symbolic-lua) in Lua.

For what it's worth, as long as he's doing his calculations in spherical coordinates, he is going to run into numerical errors at the poles -- which he probably won't overcome until he makes the switch to pseudo-Cartesian as well.  This is probably why his orbiting view is fixed to the equatorial plane.


## black hole raytrace simulator

http://christopheremoore.net/black-hole-skymap/

raytraces relativistic phenomena such as black holes and Alcubierre warp drive bubbles

does so by iterating the 4D rays of light through geodesics of arbitrary geometry.

black hole geometry is calculated via Schwarzschild metric expanded to Cartesian coordinates.
schwarzschild metric

Alcubierre warp drive bubble was already provided in Cartesian.  Way to go ADM formalism!
https://en.wikipedia.org/wiki/Alcubierre_drive

The first version of this I wrote as a luajit script.  It works with Malkia's UFO plus some
auxiliary lua libraries I have yet to post.  Find it in the lua/ folder.
That version operates via floating point FBO (much faster), yet it ran only in the view plane
and therefore had to reset all cast rays every time the camera moved.  Another downfall is you'll
find the Schwarzschild metric to be in spherical coordinate form, which caused floating point 
precision errors at the poles and along one of the dividing lines of the sphere.

The second version, and subsequently first WebGL version, is in the 'v1' folder.  It is as close
of a direct port of the Lua version and is likewise driven by FBOs, yet my current hardware
doesn't support float textures =( I'll test and rework it another time.  Am somehow withholding 
myself from writing a 24bpp fixed-point-precision raytracer.
