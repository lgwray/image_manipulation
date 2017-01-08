///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//                                              Modified:   Feng Liu
//                                              Date:       Winter 2011
//                                              Why:        Change the library file 
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <vector>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial

//computes paint_layer
void Paint_Layer(TargaImage* canvas, TargaImage* ref_image, int rad)
{
	vector<Stroke> stroke_list;
	vector<double> diff;
	int t = 25; //threshhold
	int max_x = 0;
	int max_y = 0;
	double max_diff = 0;
	double temp = 0;
	double error_area = 0;

	int rand_num = 0;

	srand(time(NULL));

	//euclidean dist for each pixel
	for (int i = 0; i < ref_image->height * ref_image->width * 4; i+=4)
	{
		diff.push_back(sqrt(pow(ref_image->data[i]-canvas->data[i], 2)
			+ pow(ref_image->data[i+1]-canvas->data[i+1], 2) 
			+ pow(ref_image->data[i+2]-canvas->data[i+2], 2)));
	}

	for (int x = 0; x < ref_image->width; x+= (2 * rad))
	{
		for (int y = 0; y < ref_image->height; y+= (2 * rad))
		{
			max_x = max_y = 0;
			error_area = max_diff = temp = 0;

			for (int i = 0; i < 2 * rad; i++)
			{
				for (int j = 0; j < 2 * rad; j++)
				{
					if ((x - rad) + j >= 0 &&
						(y - rad) + i >= 0 &&
						(x - rad) + j < ref_image->width &&
						(y - rad) + i < ref_image->height)
					{
						temp = diff[(((y - rad) + i) * ref_image->width) + ((x - rad) + j)];
						error_area += temp;

						if (temp > max_diff)
						{
							max_diff = temp;
							max_x = (x - rad) + j;
							max_y = (y - rad) + i;
						}
					}
				}
			}
			error_area /= pow((2 * rad), 2);
			//if error > t find largest error in the masked region
			// and make a stroke with that info
			if (error_area > t)
			{
				stroke_list.push_back(Stroke(rad, max_x, max_y,
					ref_image->data[((max_y * ref_image->width) + max_x) * 4], 
					ref_image->data[(((max_y * ref_image->width) + max_x) * 4) + 1],
					ref_image->data[(((max_y * ref_image->width) + max_x) * 4) + 2],
					ref_image->data[(((max_y * ref_image->width) + max_x) * 4) + 3]));
			}
		}
	}

	//paint the strokes in random order
//	random_shuffle(stroke_list.begin(), stroke_list.end());
	for (int i = (stroke_list.size() - 1); i >= 0; i--)
	{
		rand_num = rand() % (i+1);
		swap(stroke_list[i], stroke_list[rand_num]);
		canvas->Paint_Stroke(stroke_list[i]);
	}
}


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
	unsigned char temp;

	for (int i = 0; i < width * height * 4; i+=4)
	{
		temp = (data[i] * .299) + (data[i+1] * .587) + (data[i+2] * .114);
		data[i] = data[i+1] = data [i+2] = temp;
	}

//    ClearToBlack();
    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
	unsigned char red[8] = {0};
	unsigned char green[8] = {0};
	unsigned char blue[4] = {0};
	int span = 32;

	//get range of reds/greens/blues

	for(int i = 0; i < 8; i++)
	{
		if (i < 4)
		{
			red[i] += (i * span);
			green[i] += (i * span);
			blue[i] += (i * (2*span));
		}

		else
		{
			red[i] += (i * span);
			green[i] += (i * span);
		}

	}

	//map original color to the new  small collor
	for (int i = 0; i < width * height * 4; i+=4)
	{
		data[i] = red[data[i]/32];
		data[i+1] = green[data[i+1]/32];
		data[i+2] = blue[data[i+2]/64] ;
	}

//    ClearToBlack();
    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
	unsigned char red[32] = {0};
	unsigned char green[32] = {0};
	unsigned char blue[32] = {0};
	int span = 8;
	int histogram[32][32][32] = {{{0}}};
	unsigned char new_red[256] = {0};
	unsigned char new_green[256] = {0};
	unsigned char new_blue[256] = {0};

	//get range of reds/greens/blues for 32 levels for each primary

	for(int i = 0; i < 32; i++)
	{
		red[i] = (i * span);
		green[i] = (i * span);
		blue[i] = (i * span);
	}

	//quantize down to the 32 levels (for each primary)
	for (int i = 0; i < width * height * 4; i+=4)
	{
		data[i] = red[data[i]/8];
		data[i+1] = green[data[i+1]/8];
		data[i+2] = blue[data[i+2]/8] ;
	}


	//fill histogram to find most popular colors
	for (int i = 0; i < width * height * 4; i+=4)
	{
		histogram[data[i]/8][data[i+1]/8][data[i+2]/8]++;

	}

	int numbers[32*32*32] = {0};
	int numbers_pos = 0;
	for (int i = 0; i < 32; i++)
	{
		for (int j = 0; j < 32; j++)
		{
			for (int k = 0; k < 32; k++)
			{
				numbers[numbers_pos] = histogram[i][j][k];
				numbers_pos++;
			}
		}
	}	

	std::sort(std::begin(numbers), std::end(numbers));
	std::reverse(std::begin(numbers), std::end(numbers));

	//get rgb for the highest counts

	bool check = false;
	for (int x = 0; x < 256; x++)
	{
		check = false;
		for (int i = 0; i < 32; i++)
		{
			for (int j = 0; j < 32; j++)
			{
				for (int k = 0; k < 32; k++)
				{
					if(histogram[i][j][k] == numbers[x])
					{
						check = true;
						new_red[x] = red[i];
						new_green[x] = green[j];
						new_blue[x] = blue[k];
						histogram[i][j][k] = NULL;
						break;
					}
				}
				if(check)
				{
					break;
				}
			}
			if(check)
			{
				break;
			}
		
		}
	}	

	//map original color to closest chosen color

	int color = 0;
	double dist = 0;
	float check_dist = 0;
	for (int i = 0; i < height * width * 4; i+=4)
	{
		dist = sqrt(pow(data[i]-new_red[0], 2)
			+ pow(data[i+1]-new_green[0], 2) 
			+ pow(data[i+2]-new_blue[0], 2));

		color = 0;

		for (int j = 1; j < 256; j++)
		{
			check_dist = sqrt(pow(data[i]-new_red[j], 2) 
				+ pow(data[i+1]-new_green[j], 2) 
				+ pow(data[i+2]-new_blue[j], 2));

			if (check_dist < dist)
			{
				color = j;
				dist = check_dist;
			}
		}

		data[i] = new_red[color];
		data[i+1] = new_green[color];
		data[i+2] = new_blue[color];
	}


//    ClearToBlack();
    return true;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
	float* data2 = new float[width * height * 4];

	To_Grayscale();
	
	for (int i = 0; i < width * height * 4; i+=4)
	{
		data2[i] = float(data[i])/float(256);
		if(data2[i] < 0.5)
		{
			data[i] = data[i+1] = data[i+2] = 0;
		}
		else
		{
			data[i] = data[i+1] = data[i+2] = 255;
		}
	}

	delete[] data2;

    //ClearToBlack();
    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
	srand(time(0));

	float* data2 = new float[width * height * 4];
	float* order = new float[width * height];
	float avg = 0;
	float thresh = 0;
	float rand_num = 0;

	To_Grayscale();

	//get intensity array and and order array of all intensities
	for (int i = 0; i < width * height * 4; i+=4)
	{
		//add random "noise" to each pixel
		rand_num = float(rand() % 4000 - 2000) / 10000;
		data2[i] = data2[i+1] = data2[i+2] = (float(data[i])/float(256)) + rand_num;
		order[i/4] = (float(data[i])/float(256)) + rand_num;

		avg = avg + data2[i];
	}

	//calc avg intensity
	avg /= (width * height);
	
	//order the intensities in the order array
	std::sort(order, order + (width * height));

	//to preserve brightness, pick thresh where the thresh 
	//is bigger than top 1-avg% of pixels and lower than the 
	//bottom avg% of pixels
	thresh = order[int((1-avg) * width * height)];

	for (int i = 0; i < width * height * 4; i+=4)
	{
		if(data2[i] < thresh)
		{
			data[i] = data[i+1] = data[i+2] = 0;
		}
		else
		{
			data[i] = data[i+1] = data[i+2] = 255;
		}
	}

	delete[] data2;
	delete[] order;
	
    //ClearToBlack();
    return true;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
	float* data2 = new float[width * height * 4];
	float** matrix = new float*[height];
	for (int i = 0; i < height; i++) matrix[i] = new float[width];

	To_Grayscale();

	//get data in values b/w 0 and 1
	for (int i = 0; i < width * height * 4; i+=4)
	{
		data2[i] = data2[i+1] = data2[i+2] = float(data[i])/float(256);
	}

	//get data in matrix for 0-1 vals
	int pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			matrix[i][j] = data2[pos];
			pos += 4;
		}
	}
	
	//update the matrix with the FS alg
	float orig = 0;
	float error = 0;
	for (int i = 0; i < height; i++)
	{
		if (!(i % 2))
		{
			for (int j = 0; j < width; j++)
			{
				orig = matrix[i][j];
				if (orig <= .5)
				{
					matrix[i][j] = 0;
					error = orig;
				}
				else
				{
					matrix[i][j] = 255;
					error = orig - 1;
				}
				if (i + 1 < height && j - 1 > 0)
					matrix[i+1][j-1] += (error * float(float(3)/float(16)));
				if (i + 1 < height)
					matrix[i+1][j] += (error * float(float(5)/float(16)));
				if (i + 1 < height && j + 1 < width)
					matrix[i+1][j+1] += (error * float(float(1)/float(16)));
				if (j + 1 < width)
					matrix[i][j+1] += (error * float(float(7)/float(16)));
			}
		}

		else
		{
			for (int j = width - 1; j >= 0; j--)
			{
				orig = matrix[i][j];
				if (orig <= .5)
				{
					matrix[i][j] = 0; 
					error = orig;
				}
				else
				{
					matrix[i][j] = 255;
					error = orig - 1;
				}
				if (i + 1 < height && j + 1 < width)
					matrix[i+1][j+1] += (error * float(float(3)/float(16)));
				if (i + 1 < height)
					matrix[i+1][j] += (error * float(float(5)/float(16)));
				if (i + 1 < height && j - 1 > 0)
					matrix[i+1][j-1] += (error * float(float(1)/float(16)));
				if (j - 1 > 0)
					matrix[i][j-1] += (error * float(float(7)/float(16)));
			}
		}
	}

	pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			data[pos] = data[pos+1] = data[pos+2] = matrix[i][j];
			pos += 4;
		}
	}


	delete [] data2;
	for(int i = 0; i < height; i++)
		delete[] matrix[i];
	delete[] matrix;

 //   ClearToBlack();
    return true;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
	float* data2 = new float[width * height * 4];
	float* order = new float[width * height];
	float avg = 0;
	float thresh = 0;

	To_Grayscale();

	//get intensity array and and order array of all intensities
	for (int i = 0; i < width * height * 4; i+=4)
	{
		data2[i] = data2[i+1] = data2[i+2] = float(data[i])/float(256);
		order[i/4] = float(data[i])/float(256);

		avg = avg + data2[i];
	}

	//calc avg intensity
	avg /= (width * height);
	
	//order the intensities in the order array
	std::sort(order, order + (width * height));

	//to preserve brightness, pick thresh where the thresh 
	//is bigger than top 1-avg% of pixels and lower than the 
	//bottom avg% of pixels
	thresh = order[int((1-avg) * width * height)];

	for (int i = 0; i < width * height * 4; i+=4)
	{
		if(data2[i] < thresh)
		{
			data[i] = data[i+1] = data[i+2] = 0;
		}
		else
		{
			data[i] = data[i+1] = data[i+2] = 255;
		}
	}

	delete[] data2;
	delete[] order;

    //ClearToBlack();
    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
	float mask[4][4] = 
		{{0.7500, 0.3750, 0.6250, 0.2500},
		{0.0625, 1.0000, 0.8750, 0.4375},
		{0.5000, .08125, 0.9375, 0.1250},
		{0.1875, .05625, 0.3125, 0.6875}};
	float* data2 = new float[width * height * 4];
	float** matrix = new float*[height];
	for (int i = 0; i < height; i++) matrix[i] = new float[width];


	To_Grayscale();

	//get intensity array of all intensities
	for (int i = 0; i < width * height * 4; i+=4)
	{
		data2[i] = data2[i+1] = data2[i+2] = float(data[i])/float(256);
	}

	//fill matrix with intensity data to compare with mask
	int pos = 0; 
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			matrix[i][j] = data2[pos];
			pos+=4;
		}
	}

	//compare intensity of pixels with the matrix and update data
	for (int i = 0; i < height; i+=4)
	{
		for (int j = 0; j < width; j+=4)
		{
			for (int k = 0; k < 4; k++)
			{
				for (int l = 0; l < 4; l++)
				{
					if (i+k < height && j+l < width)
					{
						if (matrix[i+k][j+l] >= mask[k][l])
						{
							matrix[i+k][j+l] = 255;
						}
						else
						{
							matrix[i+k][j+l] = 0;
						}
					}
				}
			}
		}
	}

	//update data with the matrix comparisons
	pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{

			data[pos] = data[pos+1] = data[pos+2] = matrix[i][j];
			pos+=4;
//			cout<<int(data[(width * height * 4)-8])<<"\n";
//			cout<<matrix[i][j]<<"-----"<<int(data[pos])<<"pos = "<<pos<<"\n";
		}
	}
//	cout<<"pos = "<<pos;

//	for (int i = 0; i<width*height*4; i+=4) cout<<int(data[i])<<" ";
//	cout<<int(data[(width * height * 4)-8])<<"n";
	
	delete[] data2;
	for (int i = 0; i < height; i++) delete [] matrix[i];
	delete[] matrix;

//    ClearToBlack();
    return true;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
	unsigned char red[8] = {0, 36, 73, 109, 146, 182, 219, 255};
	unsigned char green[8] = {0, 36, 73, 109, 146, 182, 219, 255};
	unsigned char blue[4] = {0, 85, 170, 255};
	float** matrix = new float*[height];
	for (int i = 0; i < height; i++)  matrix[i] = new float[width * 3];


	//get data in matrix form
	int pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 3; j+=3)
		{
			matrix[i][j] = data[pos];      //red
			matrix[i][j+1] = data[pos+1];  //green
			matrix[i][j+2] = data[pos+2];  //blue
			pos += 4;                      //update data index
		}
	}

	//use fs color alg
	float orig_r, orig_g, orig_b = 0;
	float error_r, error_g, error_b = 0;
	int index_r, index_g, index_b = 0;
	float min_r, min_g, min_b = 0;

	for (int i = 0; i < height; i++)
	{
		//zig left to right on even rows
		//and zag right to left on odd rows
		if (!(i % 2))
		{
			for (int j = 0; j < width * 3; j += 3)
			{
				orig_r = matrix[i][j];
				orig_g = matrix[i][j+1];
				orig_b = matrix[i][j+2];

				if (orig_r < 0)
				{
					matrix[i][j] = red[0];
				}
				else if (orig_r > 256)
				{
					matrix[i][j] = red[7];
				}
				else
				{
					matrix[i][j] = red[int(orig_r)/32];
				}

				if (orig_g < 0)
				{
					matrix[i][j+1] = green[0];
				}
				else if (orig_g > 256)
				{
					matrix[i][j+1] = green[7];
				}
				else
				{
					matrix[i][j+1] = green[int(orig_g)/32];
				}

				if (orig_b < 0)
				{
					matrix[i][j+2] = blue[0];
				}
				else if (orig_b > 256)
				{
					matrix[i][j+2] = blue[3];
				}
				else
				{
					matrix[i][j+2] = blue[int(orig_b)/64];
				}
//
				error_r = orig_r - matrix[i][j];
				error_g = orig_g - matrix[i][j+1];
				error_b = orig_b - matrix[i][j+2];
//
				//red error distribution
				if (i + 1 < height && j - 3 > 0)
					matrix[i+1][j-3] += (error_r * float(float(3)/float(16)));
				if (i + 1 < height)
					matrix[i+1][j] += (error_r * float(float(5)/float(16)));
				if (i + 1 < height && j + 3 < width * 3)
					matrix[i+1][j+3] += (error_r * float(float(1)/float(16)));
				if (j + 3 < width *3 )
					matrix[i][j+3] += (error_r * float(float(7)/float(16)));

				//green error distribution
				if (i + 1 < height && j - 2 > 0)
					matrix[i+1][j-2] += (error_g * float(float(3)/float(16)));
				if (i + 1 < height)
					matrix[i+1][j+1] += (error_g * float(float(5)/float(16)));
				if (i + 1 < height && j + 4 < width * 3)
					matrix[i+1][j+4] += (error_g * float(float(1)/float(16)));
				if (j + 4 < width * 3)
					matrix[i][j+4] += (error_g * float(float(7)/float(16)));

				//blue error distribution
				if (i + 1 < height && j - 1 > 0)
					matrix[i+1][j-1] += (error_b * float(float(3)/float(16)));
				if (i + 1 < height)
					matrix[i+1][j+2] += (error_b * float(float(5)/float(16)));
				if (i + 1 < height && j + 5 < width * 3)
					matrix[i+1][j+5] += (error_b * float(float(1)/float(16)));
				if (j + 5 < width * 3)
					matrix[i][j+5] += (error_b * float(float(7)/float(16)));
//
			}
		}

		else
		{
			//zag right to left on odd i (height)
			for (int j = (width * 3) - 1; j >= 0; j -= 3)
			{
				orig_r = matrix[i][j-2];
				orig_g = matrix[i][j-1];
				orig_b = matrix[i][j];

				if (orig_r < 0)
				{
					matrix[i][j-2] = red[0];
				}
				else if (orig_r > 256)
				{
					matrix[i][j-2] = red[7];
				}
				else
				{
					matrix[i][j-2] = red[int(orig_r)/32];
				}

				if (orig_g < 0)
				{
					matrix[i][j-1] = green[0];
				}
				else if (orig_g > 256)
				{
					matrix[i][j-1] = green[7];
				}
				else
				{
					matrix[i][j-1] = green[int(orig_g)/32];
				}

				if (orig_b < 0)
				{
					matrix[i][j] = blue[0];
				}
				else if (orig_b > 256)
				{
					matrix[i][j] = blue[3];
				}
				else
				{
					matrix[i][j] = blue[int(orig_b)/64];
				}

//
				error_r = orig_r - matrix[i][j-2];
				error_g = orig_g - matrix[i][j-1];
				error_b = orig_b - matrix[i][j];
//
				//red error distribution
				if (i + 1 < height && j + 1 < width * 3)
					matrix[i+1][j+1] += (error_r * float(float(3)/float(16)));
				if (i + 1 < height)
					matrix[i+1][j-2] += (error_r * float(float(5)/float(16)));
				if (i + 1 < height && j - 5 > 0)
					matrix[i+1][j-5] += (error_r * float(float(1)/float(16)));
				if (j - 5 > 0)
					matrix[i][j-5] += (error_r * float(float(7)/float(16)));

				//green error distribution
				if (i + 1 < height && j + 2 < width * 3)
					matrix[i+1][j+2] += (error_g * float(float(3)/float(16)));
				if (i + 1 < height)
					matrix[i+1][j-1] += (error_g * float(float(5)/float(16)));
				if (i + 1 < height && j - 4 > 0)
					matrix[i+1][j-4] += (error_g * float(float(1)/float(16)));
				if (j - 4 > 0)
					matrix[i][j-4] += (error_g * float(float(7)/float(16)));

				//blue error distribution
				if (i + 1 < height && j + 3 < width * 3)
					matrix[i+1][j+3] += (error_b * float(float(3)/float(16)));
				if (i + 1 < height)
					matrix[i+1][j] += (error_b * float(float(5)/float(16)));
				if (i + 1 < height && j - 3 > 0)
					matrix[i+1][j-3] += (error_b * float(float(1)/float(16)));
				if (j - 3 > 0)
					matrix[i][j-3] += (error_b * float(float(7)/float(16)));
//
			}
		}
	}

	//put matrix colors into data to get new image
	pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 3; j+=3)
		{
			data[pos] = matrix[i][j];
			data[pos+1] = matrix[i][j+1];
			data[pos+2] = matrix[i][j+2];
			pos += 4;			
		}
	}

	for (int i = 0; i < height; i++)
		delete [] matrix[i];
	delete [] matrix;

//    ClearToBlack();
    return true;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

	float* a = new float[width * height * 4];
	float* b = new float[width * height * 4];
	float* c = new float[width * height * 4];

	for (int i = 0; i < width * height * 4; i+=4)
	{
		//get alpha val 0-1
		a[i+3] = float(data[i+3])/float(255);
		b[i+3] = float(pImage->data[i+3])/float(255);

		//pre-mult pixels (color * alpha)
		a[i] = data[i] * a[i+3];
		a[i+1] = data[i+1] * a[i+3];
		a[i+2] = data[i+2] * a[i+3];

		b[i] = pImage->data[i] * b[i+3];
		b[i+1] = pImage->data[i+1] * b[i+3];
		b[i+2] = pImage->data[i+2] * b[i+3];

		//calc new output image color for each pixel
		//color = color(a) = ((1 - alpha(a)) * color(b))
		c[i] = a[i] + ((1 - a[i+3]) * b[i]);  //red
		c[i+1] = a[i+1] + ((1-a[i+3]) * b[i+1]); //green
		c[i+2] = a[i+2] + ((1-a[i+3]) * b[i+2]); //blue
		c[i+3] = a[i+3] + ((1-a[i+3]) * b[i+3]); //alpha

		//divide out the alpha so not pre-mult anymore
		if (c[i+3] != 0)
		{
			c[i] = c[i]/c[i+3]; //red
			c[i+1] = c[i+1]/c[i+3]; //green
			c[i+2] = c[i+2]/c[i+3]; //blue
		}

		//get alpha as 0-255
		c[i+3] = c[i+3] * 255;

		//fill data
		data[i] = c[i];
		data[i+1] = c[i+1];
		data[i+2] = c[i+2];
		data[i+3] = c[i+3];

	}

	delete [] a;
	delete [] b;
	delete [] c;


//    ClearToBlack();
    return true;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

	float* a = new float[width * height * 4];
	float* b = new float[width * height * 4];
	float* c = new float[width * height * 4];

	for (int i = 0; i < width * height * 4; i+=4)
	{
		//get alpha val 0-1
		a[i+3] = float(data[i+3])/float(255);
		b[i+3] = float(pImage->data[i+3])/float(255);

		//pre-mult pixels (color * alpha)
		a[i] = data[i] * a[i+3];
		a[i+1] = data[i+1] * a[i+3];
		a[i+2] = data[i+2] * a[i+3];

		b[i] = pImage->data[i] * b[i+3];
		b[i+1] = pImage->data[i+1] * b[i+3];
		b[i+2] = pImage->data[i+2] * b[i+3];

		//calc new output image color for each pixel
		//color = color(a) * alpha(b)
		c[i] = a[i] * b[i+3];  //red
		c[i+1] = a[i+1] * b[i+3]; //green
		c[i+2] = a[i+2] * b[i+3]; //blue
		c[i+3] = a[i+3] * b[i+3]; //alpha

		//divide out the alpha so not pre-mult anymore
		if (c[i+3] != 0)
		{
			c[i] = c[i]/c[i+3]; //red
			c[i+1] = c[i+1]/c[i+3]; //green
			c[i+2] = c[i+2]/c[i+3]; //blue
		}

		//get alpha as 0-255
		c[i+3] = c[i+3] * 255;

		//fill data
		data[i] = c[i];
		data[i+1] = c[i+1];
		data[i+2] = c[i+2];
		data[i+3] = c[i+3];

	}

	delete [] a;
	delete [] b;
	delete [] c;

 //   ClearToBlack();
    return true;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

	float* a = new float[width * height * 4];
	float* b = new float[width * height * 4];
	float* c = new float[width * height * 4];

	for (int i = 0; i < width * height * 4; i+=4)
	{
		//get alpha val 0-1
		a[i+3] = float(data[i+3])/float(255);
		b[i+3] = float(pImage->data[i+3])/float(255);

		//pre-mult pixels (color * alpha)
		a[i] = data[i] * a[i+3];
		a[i+1] = data[i+1] * a[i+3];
		a[i+2] = data[i+2] * a[i+3];

		b[i] = pImage->data[i] * b[i+3];
		b[i+1] = pImage->data[i+1] * b[i+3];
		b[i+2] = pImage->data[i+2] * b[i+3];

		//calc new output image color for each pixel
		//color = (color(a) * (1 - alpha(b))) 
		c[i] = a[i] * (1 - b[i+3]);  //red
		c[i+1] = a[i+1] * (1 - b[i+3]); //green
		c[i+2] = a[i+2] * (1 - b[i+3]); //blue
		c[i+3] = a[i+3] * (1 - b[i+3]); //alpha

		//divide out the alpha so not pre-mult anymore
		if (c[i+3] != 0)
		{
			c[i] = c[i]/c[i+3]; //red
			c[i+1] = c[i+1]/c[i+3]; //green
			c[i+2] = c[i+2]/c[i+3]; //blue
		}

		//get alpha as 0-255
		c[i+3] = c[i+3] * 255;

		//fill data
		data[i] = c[i];
		data[i+1] = c[i+1];
		data[i+2] = c[i+2];
		data[i+3] = c[i+3];

	}

	delete [] a;
	delete [] b;
	delete [] c;

//    ClearToBlack();
    return true;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

	float* a = new float[width * height * 4];
	float* b = new float[width * height * 4];
	float* c = new float[width * height * 4];

	for (int i = 0; i < width * height * 4; i+=4)
	{
		//get alpha val 0-1
		a[i+3] = float(data[i+3])/float(255);
		b[i+3] = float(pImage->data[i+3])/float(255);

		//pre-mult pixels (color * alpha)
		a[i] = data[i] * a[i+3];
		a[i+1] = data[i+1] * a[i+3];
		a[i+2] = data[i+2] * a[i+3];

		b[i] = pImage->data[i] * b[i+3];
		b[i+1] = pImage->data[i+1] * b[i+3];
		b[i+2] = pImage->data[i+2] * b[i+3];

		//calc new output image color for each pixel
		//color = (color(a) * alpha(b)) + (color(b) * (1 - alpha(a))) 
		c[i] = (a[i] * b[i+3]) + (b[i] * (1 - a[i+3]));  //red
		c[i+1] = (a[i+1] * b[i+3]) + (b[i+1] * (1 - a[i+3])); //green
		c[i+2] = (a[i+2] * b[i+3]) + (b[i+2] * (1 - a[i+3])); //blue
		c[i+3] = (a[i+3] * b[i+3]) + (b[i+3] * (1 - a[i+3])); //alpha

		//divide out the alpha so not pre-mult anymore
		if (c[i+3] != 0)
		{
			c[i] = c[i]/c[i+3]; //red
			c[i+1] = c[i+1]/c[i+3]; //green
			c[i+2] = c[i+2]/c[i+3]; //blue
		}

		//get alpha as 0-255
		c[i+3] = c[i+3] * 255;

		//fill data
		data[i] = c[i];
		data[i+1] = c[i+1];
		data[i+2] = c[i+2];
		data[i+3] = c[i+3];

	}

	delete [] a;
	delete [] b;
	delete [] c;

//    ClearToBlack();
    return true;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

	float* a = new float[width * height * 4];
	float* b = new float[width * height * 4];
	float* c = new float[width * height * 4];

	for (int i = 0; i < width * height * 4; i+=4)
	{
		//get alpha val 0-1
		a[i+3] = float(data[i+3])/float(255);
		b[i+3] = float(pImage->data[i+3])/float(255);

		//pre-mult pixels (color * alpha)
		a[i] = data[i] * a[i+3];
		a[i+1] = data[i+1] * a[i+3];
		a[i+2] = data[i+2] * a[i+3];

		b[i] = pImage->data[i] * b[i+3];
		b[i+1] = pImage->data[i+1] * b[i+3];
		b[i+2] = pImage->data[i+2] * b[i+3];

		//calc new output image color for each pixel
		//color = (color(a) * (1 - alpha(b))) + (color(b) * (1 - alpha(a))) 
		c[i] = (a[i] * (1 - b[i+3])) + (b[i] * (1 - a[i+3]));  //red
		c[i+1] = (a[i+1] * (1 - b[i+3])) + (b[i+1] * (1 - a[i+3])); //green
		c[i+2] = (a[i+2] * (1 - b[i+3])) + (b[i+2] * (1 - a[i+3])); //blue
		c[i+3] = (a[i+3] * (1 - b[i+3])) + (b[i+3] * (1 - a[i+3])); //alpha

		//divide out the alpha so not pre-mult anymore
		if (c[i+3] != 0)
		{
			c[i] = c[i]/c[i+3]; //red
			c[i+1] = c[i+1]/c[i+3]; //green
			c[i+2] = c[i+2]/c[i+3]; //blue
		}

		//get alpha as 0-255
		c[i+3] = c[i+3] * 255;

		//fill data
		data[i] = c[i];
		data[i+1] = c[i+1];
		data[i+2] = c[i+2];
		data[i+3] = c[i+3];

	}

	delete [] a;
	delete [] b;
	delete [] c;

//    ClearToBlack();
    return true;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
	float box[5][5];
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			box[i][j] = float(1)/float(25);
		}
	}

	float** data2 = new float*[height];
	float** temp = new float*[height];
	for (int i = 0; i < height; i++)
	{
		data2[i] = new float[width * 3];
		temp[i] = new float[width * 3];
	}

	int pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 3; j+=3)
		{
			data2[i][j] = data[pos];
			data2[i][j+1] = data[pos+1];
			data2[i][j+2] = data[pos+2];
			
			temp[i][j] = 0;
			temp[i][j+1] = 0;
			temp[i][j+2] = 0;

			pos+=4;
		}
	}



	for (int i = 2; i < height - 2; i++)
	{
		for (int j = 6; j < (width * 3) - 6; j+=3)
		{
			for (int k = 0; k < 5; k++)
			{
				for (int l = 0; l < 5; l++)
				{
					temp[i][j] += data2[k + (i-2)][(l*3) + (j-6)] * box[k][l];
					temp[i][j+1] += data2[k + (i-2)][(l*3) + (j-6) + 1] * box[k][l];
					temp[i][j+2] += data2[k + (i-2)][(l*3) + (j-6) + 2] * box[k][l];
				}
			}
		}
	}

	pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 3; j+=3)
		{
			data[pos] = temp[i][j];
			data[pos+1] = temp[i][j+1];
			data[pos+2] = temp[i][j+2];
			pos += 4;
		}
	}

	for (int i = 0; i < height; i++)
	{
		delete [] temp[i];
		delete [] data2[i];
	}
	delete[] temp;
	delete[] data2;

 //   ClearToBlack();
    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
	float bart[5][5] = {{1, 2, 3, 2, 1},
	                    {2, 4, 6, 4, 2},
	                    {3, 6, 9, 6, 3},
	                    {2, 4, 6, 4, 2},
	                    {1, 2, 3, 2, 1}};

	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			bart[i][j] = float(bart[i][j])/float(81);
		}
	}

	float** data2 = new float*[height];
	float** temp = new float*[height];
	for (int i = 0; i < height; i++)
	{
		data2[i] = new float[width * 3];
		temp[i] = new float[width * 3];
	}

	int pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 3; j+=3)
		{
			data2[i][j] = data[pos];
			data2[i][j+1] = data[pos+1];
			data2[i][j+2] = data[pos+2];
			
			temp[i][j] = 0;
			temp[i][j+1] = 0;
			temp[i][j+2] = 0;

			pos+=4;
		}
	}



	for (int i = 2; i < height - 2; i++)
	{
		for (int j = 6; j < (width * 3) - 6; j+=3)
		{
			for (int k = 0; k < 5; k++)
			{
				for (int l = 0; l < 5; l++)
				{
					temp[i][j] += data2[k + (i-2)][(l*3) + (j-6)] * bart[k][l];
					temp[i][j+1] += data2[k + (i-2)][(l*3) + (j-6) + 1] * bart[k][l];
					temp[i][j+2] += data2[k + (i-2)][(l*3) + (j-6) + 2] * bart[k][l];
				}
			}
		}
	}

	pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 3; j+=3)
		{
			data[pos] = temp[i][j];
			data[pos+1] = temp[i][j+1];
			data[pos+2] = temp[i][j+2];
			pos += 4;
		}
	}

	for (int i = 0; i < height; i++)
	{
		delete [] temp[i];
		delete [] data2[i];
	}
	delete[] temp;
	delete[] data2;

//    ClearToBlack();
    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
	float gaus[5][5] = {{1, 4, 6, 4, 1},
	                    {4, 16, 24, 16, 4},
	                    {6, 24, 36, 24, 6},
	                    {4, 16, 24, 16, 4},
	                    {1, 4, 6, 4, 1}};

	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			gaus[i][j] = float(gaus[i][j])/float(256);
		}
	}

	float** data2 = new float*[height];
	float** temp = new float*[height];
	for (int i = 0; i < height; i++)
	{
		data2[i] = new float[width * 3];
		temp[i] = new float[width * 3];
	}

	int pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 3; j+=3)
		{
			data2[i][j] = data[pos];
			data2[i][j+1] = data[pos+1];
			data2[i][j+2] = data[pos+2];
			
			temp[i][j] = 0;
			temp[i][j+1] = 0;
			temp[i][j+2] = 0;

			pos+=4;
		}
	}



	for (int i = 2; i < height - 2; i++)
	{
		for (int j = 6; j < (width * 3) - 6; j+=3)
		{
			for (int k = 0; k < 5; k++)
			{
				for (int l = 0; l < 5; l++)
				{
					temp[i][j] += data2[k + (i-2)][(l*3) + (j-6)] * gaus[k][l];
					temp[i][j+1] += data2[k + (i-2)][(l*3) + (j-6) + 1] * gaus[k][l];
					temp[i][j+2] += data2[k + (i-2)][(l*3) + (j-6) + 2] * gaus[k][l];
				}
			}
		}
	}

	pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 3; j+=3)
		{
			data[pos] = temp[i][j];
			data[pos+1] = temp[i][j+1];
			data[pos+2] = temp[i][j+2];
			pos += 4;
		}
	}

	for (int i = 0; i < height; i++)
	{
		delete [] temp[i];
		delete [] data2[i];
	}
	delete[] temp;
	delete[] data2;

 //   ClearToBlack();
    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
	double** gaus = new double*[N];
	for (int i = 0; i < N; i++)
		gaus[i] = new double[N];

	for (int i = 0; i < N; i++)
		gaus[0][i] = Binomial(N-1, i);

	double total = 0;
	int border = N/2;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			gaus[i][j] = gaus[0][i] * gaus[0][j];
			total += gaus[i][j];
		}
	}

	float** data2 = new float*[height];
	float** temp = new float*[height];
	for (int i = 0; i < height; i++)
	{
		data2[i] = new float[width * 3];
		temp[i] = new float[width * 3];
	}

	int pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 3; j+=3)
		{
			data2[i][j] = data[pos];
			data2[i][j+1] = data[pos+1];
			data2[i][j+2] = data[pos+2];
			
			temp[i][j] = 0;
			temp[i][j+1] = 0;
			temp[i][j+2] = 0;

			pos+=4;
		}
	}
	

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < (width * 3); j+=3)
		{
			for (int k = 0; k < N; k++)
			{
				for (int l = 0; l < N; l++)
				{
					if((i - border) + k >= 0 &&
						(j - (border * 3)) + (l * 3) >= 0 &&
						(i - border) + k < height &&
						(j - (border * 3)) + (l * 3) < width * 3)
					{
					temp[i][j] += data2[(i - border) + k][(j - (border * 3)) + (l * 3)] * (gaus[k][l]/total);
					temp[i][j+1] += data2[(i - border) + k][(j - (border * 3)) + (l * 3) +1] * (gaus[k][l]/total);
					temp[i][j+2] += data2[(i - border) + k][(j - (border * 3)) + (l * 3) +2] * (gaus[k][l]/total);
					}
				}
			}
		}
	}

	pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < (width * 3); j+=3)
		{
//			if (i >= border && j >= border * 3 && i < height - border && j < (width * 3) - (border * 3))
//			{
			data[pos] = temp[i][j];
			data[pos+1] = temp[i][j+1];
			data[pos+2] = temp[i][j+2];
//			}
			pos += 4;
		}
	}

	for (int i = 0; i < height; i++)
	{
		delete [] temp[i];
		delete [] data2[i];
	}
	delete[] temp;
	delete[] data2;

	for (int i = 0; i < N; i++)
		delete [] gaus[i];
	delete [] gaus;

//    ClearToBlack();
   return true;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
	float gaus[5][5] = {{1, 4, 6, 4, 1},
	                    {4, 16, 24, 16, 4},
	                    {6, 24, 36, 24, 6},
	                    {4, 16, 24, 16, 4},
	                    {1, 4, 6, 4, 1}};

	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			gaus[i][j] = float(gaus[i][j])/float(256);
		}
	}

	float** data2 = new float*[height];
	float** temp = new float*[height];
	for (int i = 0; i < height; i++)
	{
		data2[i] = new float[width * 3];
		temp[i] = new float[width * 3];
	}

	int pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 3; j+=3)
		{
			data2[i][j] = data[pos];
			data2[i][j+1] = data[pos+1];
			data2[i][j+2] = data[pos+2];
			
			temp[i][j] = 0;
			temp[i][j+1] = 0;
			temp[i][j+2] = 0;

			pos+=4;
		}
	}



	for (int i = 2; i < height - 2; i++)
	{
		for (int j = 6; j < (width * 3) - 6; j+=3)
		{
			for (int k = 0; k < 5; k++)
			{
				for (int l = 0; l < 5; l++)
				{
					temp[i][j] += data2[k + (i-2)][(l*3) + (j-6)] * gaus[k][l];
					temp[i][j+1] += data2[k + (i-2)][(l*3) + (j-6) + 1] * gaus[k][l];
					temp[i][j+2] += data2[k + (i-2)][(l*3) + (j-6) + 2] * gaus[k][l];
				}
			}
		}
	}

	pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 3; j+=3)
		{
			if(float(data[pos]) - temp[i][j]  < 0) data[pos] = 0;
			else if(float(data[pos]) - temp[i][j] > 255) data[pos] = 255;
			else data[pos] -= temp[i][j];
			
			if(float(data[pos+1]) - temp[i][j+1]  < 0) data[pos+1] = 0;
			else if(float(data[pos+1]) - temp[i][j] > 255) data[pos+1] = 255;
			else data[pos+1] -= temp[i][j+1];
			
			if(float(data[pos+2]) - temp[i][j+2]  < 0) data[pos+2] = 0;
			else if(float(data[pos+2]) - temp[i][j] > 255) data[pos+2] = 255;
			else data[pos+2] -= temp[i][j+2];
			
			pos += 4;
		}
	}

	for (int i = 0; i < height; i++)
	{
		delete [] temp[i];
		delete [] data2[i];
	}
	delete[] temp;
	delete[] data2;

    //ClearToBlack();
    return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
	float gaus[5][5] = {{1, 4, 6, 4, 1},
	                    {4, 16, 24, 16, 4},
	                    {6, 24, 36, 24, 6},
	                    {4, 16, 24, 16, 4},
	                    {1, 4, 6, 4, 1}};

	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			gaus[i][j] = float(gaus[i][j])/float(256);
		}
	}

	float** data2 = new float*[height];
	float** temp = new float*[height];
	for (int i = 0; i < height; i++)
	{
		data2[i] = new float[width * 3];
		temp[i] = new float[width * 3];
	}

	int pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 3; j+=3)
		{
			data2[i][j] = data[pos];
			data2[i][j+1] = data[pos+1];
			data2[i][j+2] = data[pos+2];
			
			temp[i][j] = 0;
			temp[i][j+1] = 0;
			temp[i][j+2] = 0;

			pos+=4;
		}
	}



	for (int i = 2; i < height - 2; i++)
	{
		for (int j = 6; j < (width * 3) - 6; j+=3)
		{
			for (int k = 0; k < 5; k++)
			{
				for (int l = 0; l < 5; l++)
				{
					temp[i][j] += data2[k + (i-2)][(l*3) + (j-6)] * gaus[k][l];
					temp[i][j+1] += data2[k + (i-2)][(l*3) + (j-6) + 1] * gaus[k][l];
					temp[i][j+2] += data2[k + (i-2)][(l*3) + (j-6) + 2] * gaus[k][l];
				}
			}
		}
	}

	pos = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 3; j+=3)
		{
			if(float(data[pos]) + (float(data[pos]) - temp[i][j]) < 0) 
				data[pos] = 0;
			else if(float(data[pos]) + (float(data[pos]) - temp[i][j]) > 255) 
				data[pos] = 255;
			else 
				data[pos] = data[pos] + (data[pos] - temp[i][j]);
			
			if(float(data[pos+1]) + (float(data[pos+1]) - temp[i][j+1]) < 0) 
				data[pos+1] = 0;
			else if(float(data[pos+1]) + (float(data[pos+1]) - temp[i][j+1]) > 255)
				data[pos+1] = 255;
			else 
				data[pos+1] = data[pos+1] + (data[pos+1] - temp[i][j+1]);
			
			if(float(data[pos+2]) + (float(data[pos+2]) - temp[i][j+2]) < 0) 
				data[pos+2] = 0;
			else if(float(data[pos+2]) + (float(data[pos+2]) - temp[i][j+2]) > 255) 
				data[pos+2] = 255;
			else data[pos+2] = data[pos+2] + (data[pos+2] - temp[i][j+2]);
			
			pos += 4;
		}
	}

	for (int i = 0; i < height; i++)
	{
		delete [] temp[i];
		delete [] data2[i];
	}
	delete[] temp;
	delete[] data2;

 //   ClearToBlack();
    return true;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
	TargaImage* canvas = new TargaImage(width, height);
	TargaImage* ref_image;

	for (int i = 0; i < width * height * 4; i+=4)
		canvas->data[i] = canvas->data[i+1] = canvas->data[i+2] = canvas->data[i+3] = 0;

	int rad[3] = {7, 3, 1};

	for (int i = 0; i < 3; i++)
	{
		//get ref image as gauss filter size R
		ref_image = new TargaImage(width, height, data);
		ref_image->Filter_Gaussian_N((2 * rad[i]) + 1);

		//paint a layer
		Paint_Layer(canvas, ref_image, rad[i]);

	}

	for (int i = 0; i < width * height * 4; i++)
		data[i] = canvas->data[i];


//    ClearToBlack();
    return true;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
	int new_height = height/2;
	int new_width = width/2;
	int pos = 0;

	float** data2 = new float*[height];
	for (int i = 0; i < height; i++) data2[i] = new float[width * 4];

	float ** data_small = new float*[new_height];
	for (int i = 0; i < new_height; i++) data_small[i] = new float[new_width * 4];

	float mask[3][3] = {{float(1)/float(16), float(1)/float(8), float(1)/float(16)},
						{float(1)/float(8), float(1)/float(4), float(1)/float(8)},
						{float(1)/float(16), float(1)/float(8), float(1)/float(16)}};


	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 4; j+=4)
		{
			data2[i][j] = data[pos];
			data2[i][j+1] = data[pos+1];
			data2[i][j+2] = data[pos+2];
			data2[i][j+3] = data[pos+3];
			pos += 4;
		}
	}


	for (int i = 0; i < new_height; i++)
	{
		for (int j = 0; j < new_width * 4; j+=4)
		{
			data_small[i][j] = data_small[i][j+1] = data_small[i][j+2] = data_small[i][j+3] = 0;
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					if (((2*i) - 1) + k >= 0 &&
						((2*j) - 4) + (l*4) >= 0 &&
						((2*i) - 1) + k < height &&
						((2*j) - 4) + (l*4) < width * 4)
					{
						data_small[i][j] += 
							data2[((2*i) - 1) + k][((2*j) - 4) + (l*4)] * mask[k][l];
						data_small[i][j+1] += 
							data2[((2*i) - 1) + k][(((2*j)+1) - 4) + (l*4)] * mask[k][l];
						data_small[i][j+2] += 
							data2[((2*i) - 1) + k][(((2*j)+2) - 4) + (l*4)] * mask[k][l];
						data_small[i][j+3] += 
							data2[((2*i) - 1) + k][(((2*j)+3) - 4) + (l*4)] * mask[k][l];
						
					}
				}
			}
		}
	}

	delete [] data;
	data = new unsigned char[new_width * new_height * 4];
	height = new_height;
	width = new_width;
	pos = 0;

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 4; j+=4)
		{

			data[pos] = data_small[i][j];
			data[pos+1] = data_small[i][j+1];
			data[pos+2] = data_small[i][j+2];
			data[pos+3] = data_small[i][j+3];

			pos += 4;
		}
	}

//    ClearToBlack();
    return true;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
	float ee[3][3] = {{float(1)/float(16), float(1)/float(8), float(1)/float(16)},
						{float(1)/float(8), float(1)/float(4), float(1)/float(8)},
						{float(1)/float(16), float(1)/float(8), float(1)/float(16)}};

	float oo[4][4] = {{float(1)/float(64), float(3)/float(64), float(3)/float(64), float(1)/float(64)},
						{float(3)/float(64), float(9)/float(64), float(9)/float(64), float(3)/float(64)},
						{float(3)/float(64), float(9)/float(64), float(9)/float(64), float(3)/float(64)},
						{float(1)/float(64), float(3)/float(64), float(3)/float(64), float(1)/float(64)}};

	float eo[3][4] = {{float(1)/float(32), float(3)/float(32), float(3)/float(32), float(1)/float(32)},
						{float(2)/float(32), float(6)/float(32), float(6)/float(32), float(2)/float(32)},
						{float(1)/float(32), float(3)/float(32), float(3)/float(32), float(1)/float(32)}};

	float oe[4][3] = {{float(1)/float(32), float(2)/float(32), float(1)/float(32)},
						{float(3)/float(32), float(6)/float(32), float(3)/float(32)},
						{float(3)/float(32), float(6)/float(32), float(3)/float(32)},
						{float(1)/float(32), float(2)/float(32), float(1)/float(32)}};

	int new_height = height * 2;
	int new_width = width * 2;
	int pos = 0;

	float** new_data = new float*[new_height];
	for (int i = 0; i < new_height; i++) new_data[i] = new float[new_width * 4];

	float** data2 = new float*[height];
	for (int i = 0; i < height; i++) data2[i] = new float[width * 4];

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 4; j+=4)
		{
			data2[i][j] = data[pos];
			data2[i][j+1] = data[pos+1];
			data2[i][j+2] = data[pos+2];
			data2[i][j+3] = data[pos+3];
			pos += 4;
		}
	}

	for (int i = 0; i < new_height; i++)
	{
		for (int j = 0; j < new_width * 4; j+=4)
		{
			new_data[i][j] = new_data[i][j+1] = new_data[i][j+2] = new_data[i][j+3] = 0;

			//even i and even j
			if (!(i % 2) && !((j/4) % 2))
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						if (((i/2) - 1) + k >= 0 &&
							((j/2) - 4) + (l*4) >= 0 &&
							((i/2) - 1) + k < height &&
							((j/2) - 4) + (l*4) < width * 4)
						{
							new_data[i][j] += 
								data2[((i/2) - 1) + k][((((j/4)/2) * 4) - 4) + (l*4)] * ee[k][l];
							new_data[i][j+1] += 
								data2[((i/2) - 1) + k][(((((j/4)/2) * 4) +1) - 4) + (l*4)] * ee[k][l];
							new_data[i][j+2] += 
								data2[((i/2) - 1) + k][(((((j/4)/2) * 4) +2) - 4) + (l*4)] * ee[k][l];
							new_data[i][j+3] += 
								data2[((i/2) - 1) + k][(((((j/4)/2) * 4) +3) - 4) + (l*4)] * ee[k][l];	
						}
					}
				}
			}

			//if i odd and j odd
			if ((i % 2) && ((j/4) % 2))
			{
				for (int k = 0; k < 4; k++)
				{
					for (int l = 0; l < 4; l++)
					{
						if (((i/2) - 1) + k >= 0 &&
							((((j/4)/2) * 4) - 4) + (l*4) >= 0 &&
							((i/2) - 1) + k < height &&
							((((j/4)/2) * 4) - 4) + (l*4) < width * 4)
						{
							new_data[i][j] += 
								data2[((i/2) - 1) + k][((((j/4)/2) * 4) - 4) + (l*4)] * oo[k][l];
							new_data[i][j+1] += 
								data2[((i/2) - 1) + k][(((((j/4)/2) * 4) +1) - 4) + (l*4)] * oo[k][l];
							new_data[i][j+2] += 
								data2[((i/2) - 1) + k][(((((j/4)/2) * 4) +2) - 4) + (l*4)] * oo[k][l];
							new_data[i][j+3] += 
								data2[((i/2) - 1) + k][(((((j/4)/2) * 4) +3) - 4) + (l*4)] * oo[k][l];	
						}
					}
				}
			}

			//if even i and odd j
			if (!(i % 2) && ((j/4) % 2))
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 4; l++)
					{
						if (((i/2) - 1) + k >= 0 &&
							((j/2) - 4) + (l*4) >= 0 &&
							((i/2) - 1) + k < height &&
							((j/2) - 4) + (l*4) < width * 4)
						{
							new_data[i][j] += 
								data2[((i/2) - 1) + k][((((j/4)/2) * 4) - 4) + (l*4)] * eo[k][l];
							new_data[i][j+1] += 
								data2[((i/2) - 1) + k][(((((j/4)/2) * 4) +1) - 4) + (l*4)] * eo[k][l];
							new_data[i][j+2] += 
								data2[((i/2) - 1) + k][(((((j/4)/2) * 4) +2) - 4) + (l*4)] * eo[k][l];
							new_data[i][j+3] += 
								data2[((i/2) - 1) + k][(((((j/4)/2) * 4) +3) - 4) + (l*4)] * eo[k][l];	
						}
					}
				}
			}

			//if odd i even j
			if ((i % 2) && !((j/4) % 2))
			{
				for (int k = 0; k < 4; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						if (((i/2) - 1) + k >= 0 &&
							((j/2) - 4) + (l*4) >= 0 &&
							((i/2) - 1) + k < height &&
							((j/2) - 4) + (l*4) < width * 4)
						{
							new_data[i][j] += 
								data2[((i/2) - 1) + k][((((j/4)/2) * 4) - 4) + (l*4)] * oe[k][l];
							new_data[i][j+1] += 
								data2[((i/2) - 1) + k][(((((j/4)/2) * 4) +1) - 4) + (l*4)] * oe[k][l];
							new_data[i][j+2] += 
								data2[((i/2) - 1) + k][(((((j/4)/2) * 4) +2) - 4) + (l*4)] * oe[k][l];
							new_data[i][j+3] += 
								data2[((i/2) - 1) + k][(((((j/4)/2) * 4) +3) - 4) + (l*4)] * oe[k][l];	
						}
					}
				}
			}
		}
	}


	for (int i = 0; i < height; i++)
		delete [] data2[i];
	delete [] data2;


	delete [] data;
	data = new unsigned char[new_width * new_height * 4];
	height = new_height;
	width = new_width;
	pos = 0;

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width * 4; j+=4)
		{
			data[pos] = new_data[i][j];
			data[pos+1] = new_data[i][j+1];
			data[pos+2] = new_data[i][j+2];
			data[pos+3] = new_data[i][j+3];
			pos += 4;

		}
	}


	for (int i = 0; i < height; i++)
		delete [] new_data[i];
	delete [] new_data;



//    ClearToBlack();
    return true;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
	

    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
	//float mask[4][4] = {{float(1)/float(64), float(3)/float(64), float(3)/float(64), float(1)/float(64)},
	//					{float(3)/float(64), float(9)/float(64), float(9)/float(64), float(3)/float(64)},
	//					{float(3)/float(64), float(9)/float(64), float(9)/float(64), float(3)/float(64)},
	//					{float(1)/float(64), float(3)/float(64), float(3)/float(64), float(1)/float(64)}};

	//int pos = 0;
	//int orig_x, orig_y = 0;
	//float pi = 3.14159265;

	//float** data_r = new float*[height];
	//for (int i = 0; i < height; i++) data_r[i] = new float[width];

	//float** data_g = new float*[height];
	//for (int i = 0; i < height; i++) data_g[i] = new float[width];

	//float** data_b = new float*[height];
	//for (int i = 0; i < height; i++) data_b[i] = new float[width];

	//float** data_a = new float*[height];
	//for (int i = 0; i < height; i++) data_a[i] = new float[width];

	//float** new_data_r = new float*[height];
	//for (int i = 0; i < height; i++) new_data_r[i] = new float[width];

	//float** new_data_g = new float*[height];
	//for (int i = 0; i < height; i++) new_data_g[i] = new float[width];

	//float** new_data_b = new float*[height];
	//for (int i = 0; i < height; i++) new_data_b[i] = new float[width];

	//float** new_data_a = new float*[height];
	//for (int i = 0; i < height; i++) new_data_a[i] = new float[width];

	//for (int i = 0; i < height; i++)
	//{
	//	for (int j = 0; j < width; j++)
	//	{
	//		data_r[i][j] = data[pos];
	//		data_g[i][j] = data[pos+1];
	//		data_b[i][j] = data[pos+2];
	//		data_a[i][j] = data[pos+3];
	//		pos += 4;
	//	}
	//}


	//for (int i = 0; i < height; i++)
	//{
	//	for (int j = 0; j < width; j++)
	//	{
	//		new_data_r[i][j] = new_data_g[i][j] = new_data_b[i][j] = new_data_a[i][j] = 0;
	//		for (int k = 0; k < 3; k++)
	//		{
	//			for (int l = 0; l < 3; l++)
	//			{
	//				/*orig_y = ((sin((angleDegrees * pi)/180)*i) + (cos((angleDegrees * pi)/180)*j));
	//				orig_x = ((cos((angleDegrees * pi)/180)*j) - (sin((angleDegrees * pi)/180)*i));*/
	//				orig_y = ((sin(angleDegrees)*i) + (cos(angleDegrees)*j));
	//				orig_x = ((cos(angleDegrees)*j) - (sin(angleDegrees)*i));
	//				if ((orig_y - 1) + k >= 0 &&
	//					(orig_x - 1) + l >= 0 &&
	//					(orig_y - 1) + k < height &&
	//					(orig_x - 1) + l < width)
	//				{
	//					new_data_r[i][j] +=
	//						data_r[(orig_y - 1) + k][(orig_x - 1) + l] * mask[k][j];
	//					new_data_g[i][j] +=
	//						data_g[(orig_y - 1) + k][(orig_x - 1) + l] * mask[k][j];
	//					new_data_b[i][j] +=
	//						data_b[(orig_y - 1) + k][(orig_x - 1) + l] * mask[k][j];
	//					new_data_a[i][j] +=
	//						data_a[(orig_y - 1) + k][(orig_x - 1) + l] * mask[k][j];
	//				}
	//			}
	//		}
	//	}
	//}

	//pos = 0;
	//for (int i = 0; i < height; i++)
	//{
	//	for (int j = 0; j < width; j++)
	//	{
	//		data[pos] = new_data_r[i][j];
	//		data[pos+1] = new_data_g[i][j];
	//		data[pos+2] = new_data_b[i][j];
	//		data[pos+3] = new_data_a[i][j];
	//		pos += 4;

	//	}
	//}


	//for (int i = 0; i < height; i++)
	//{
	//	delete [] new_data_r[i];
	//	delete [] new_data_g[i];
	//	delete [] new_data_b[i];
	//	delete [] new_data_a[i];
	//	delete [] data_r[i];
	//	delete [] data_g[i];
	//	delete [] data_b[i];
	//	delete [] data_a[i];
	//}
	//delete [] new_data_r;
	//delete [] new_data_g;
	//delete [] new_data_b;
	//delete [] new_data_a;
	//delete [] data_r;
	//delete [] data_g;
	//delete [] data_b;
	//delete [] data_a;
	
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}


