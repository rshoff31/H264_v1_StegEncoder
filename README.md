# H264_v1_StegEncoder
Program used in Audio/Video Compression Algorithms class.

########## - H264_v1 - Steganography Edition - ##########
Use the following code as main: A00_Main.m

## NOTES ##
- This program does a JSTEG style Steganography embedding of text file data into coefficients with a value of -2, -1, 2, and 3:
   -If message bit = 0, and coefficient = -2, then no change.
   -If message bit = 1, and coefficient = -2, then change coefficient to -1.
   -If message bit = 0, and coefficient = -1, then change coefficient to -2.
   -If message bit = 1, and coefficient = -1, then no change.
   -If message bit = 0, and coefficient = 2, then no change.
   -If message bit = 1, and coefficient = 2, then change coefficient to 3.
   -If message bit = 0, and coefficient = 3, then change coefficient to 2.
   -If message bit = 1, and coefficient = 3, then no change.

- This program does not do decoding of data; the purpose of the program is simply to show how the statistics of embedding during the encoding process affect the statistics of the image.

## Usage ##
1. Open A00_Main.m with Matlab.
2. Modify variables as desired under 'Variables to Modify' section.
3. Run Code!!!
4. OUTPUT
    a.) Encoder:
        i.)   Large window, showing original, color video.
        ii.)  Small window, showing Grayscale video that is compressed.
        iii.) In console, first diff block, DCT block, Quantize DCT block, and Quantized DCT block modified with secret message shown, and H.264 video statistics.
    b.) Decoder:
        i.)  Large window, showing the following videos:
              Top: original video frames.
              Mid: Reconstructed, quantized DCT frames, without modification.
              Bot: Reconstructed, quantized DCT frames, with modification.
        ii.) In console, show reconstructed 
