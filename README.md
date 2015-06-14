Due to our inability to make video works in eclipse, we made a version wihthout
video.
There is a method in TangibleGame (changeInputImage) which should allow to use
our code with a video if we pass in arguments the image of each frame.

The class ImageProcessing is not really used, on some machine we can use it to  
read the video with the hough and the sobel. Not sure it works on linux or 
macosx

We decided to put in folder cs211/ not only the .java which are in package game
but also the source images and other things we needed.

in folder TGame, we have the full eclipse project with all the libraries needed
but there is a few absolute paths who sould need to be changed.