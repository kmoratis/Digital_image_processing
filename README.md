# Digital_image_processing

#### Raw to RGB image transform

 - Reading Raw image (.DNG) and extracting useful metadata

 - White-balancing, adjusted to Bayer pattern ( 'RGGB', 'BGGR', 'GRBG', 'GBRG' )

 - Nearest-neighbor interpolation - bilinear interpolation

 - Space transforms ( Cam, XYZ, sRGB )

#### Optical Character Recognition

- Find image rotation angle of a text image, using its DFT and horizontal projection

- Image rotation implementation and white-padding

- Find the contours of a character contained in an image

- Create training dataset from given image

- Training - validation set split

- Train K-NN classification models

- Evaluate models using given test image
