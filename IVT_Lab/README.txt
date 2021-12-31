1) All cpp codes except for Session1.cpp are especified by part of the session
	Session2p4 correspond to Session 2 item 4
	Session2p5 correspond to Session 2 item 5
	.............
	.............
2) Input files
All codes needs to read the file "lena_256x256.raw"
Session2p6 needs to read the file Image2lenaNoisy.raw which is the result after execute Session2p5

3) Output files
Session1-->ImageS1.raw (256x256)(32bits)

Session2p4-->Uniform.raw (256x256)(32bits), Gaussian.raw(256x256)(32bits)
Session2p5-->Image2lenaNoisy.raw (256x256)(32bits)
Session2p6-->Image2Blur.raw (256x256)(32bits)
Note: This code needs the previous image to run Image2lenaNoisy.raw

Session3p7-->identity.raw (256x256)(32bits), matrix.raw(256x256)(32bits), matrixTransp.raw(256x256)(32bits)
Session3p8-->IDCT.raw (256x256)(32bits), LenaDCT2.raw(256x256)(32bits)

Session4p9-->lena8bpp.raw(256x256)(8bits), lenajpeg.raw (256x256)(32bits)
To read lena8bpp.raw change the pixel to 8bits. 

Session5p10-->file.txt, ImageDCTerms.raw (32x32)(32bits), ImageRecostructed32x32.raw(32x32)(32bits)
file.txt contains all DC terms differences
Session5p11-->ImageRecostructedRLE (256x256)(32bits)
Session5p12-->ImmageRecostructedRLE (256x256)(32bits)