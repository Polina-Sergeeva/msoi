#include <vector>
#include <windows.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>

int height = 0, width = 0, height_temp = 0, width_temp = 0;

std::vector<float**> get_rgb_components(FILE* in);
float correlation(float** A, float** B);
void autocorrelation(float** component, std::string component_name);
float** get_sample(float** component, int y, int x);
std::vector<float**> rgb_to_ycbcr(std::vector<float**> rgb, FILE* in);
std::vector<float**> ycbcr_to_rgb(std::vector<float**> ycbcr);
float PSNR(float** original, float** recovered);
float** decimation1(float** component, int coef);
float** decimation2(float** component, int coef);
float** decimation_recover(float** component, int coef);
void histogram(float** component, std::string component_name, bool is_diff);
void enthropy(std::map<float, unsigned int>& values, std::string& component_name);
float** differential_coding(float** component);
std::vector<float**> get_subframes(float** component, FILE* in);
void subframes_modification(std::vector<float**>& subframes);

int main(int argc, char const *argv[])
{
	FILE* in = fopen(argv[argc - 1], "rb");

	//2-3
	std::vector<float**> components_rgb = get_rgb_components(in);	

	//4.a
	printf("\ncorrelations coeffitients:\nr[R,G]=%f\nr[R,B]=%f\nr[G,B]=%f\n", correlation(components_rgb[0], components_rgb[1]), 
						  correlation(components_rgb[0], components_rgb[2]), correlation(components_rgb[1], components_rgb[2]));
	//4.b
	autocorrelation(components_rgb[0], "R");
	autocorrelation(components_rgb[1], "G");
	autocorrelation(components_rgb[2], "B");
	printf("\nautocorrelation functions for rgb done!\n");

	//6
	std::vector<float**> components_ycbcr = rgb_to_ycbcr(components_rgb, in);

	//5(4.a)
	printf("\ncorrelations coeffitients:\nr[Y,Cb]=%f\nr[Y,Cr]=%f\nr[Cb,Cr]=%f\n",correlation(components_ycbcr[0], components_ycbcr[1]), 
						 correlation(components_ycbcr[0], components_ycbcr[2]), correlation(components_ycbcr[1], components_ycbcr[2]));
	//5(4.b)
	autocorrelation(components_ycbcr[0], "Y");
	autocorrelation(components_ycbcr[1], "Cb");
	autocorrelation(components_ycbcr[2], "Cr");
	printf("\nautocorrelation functions for ycbcr done!\n");

	//7
	std::vector<float**> components_rgb_recovered = ycbcr_to_rgb(components_ycbcr);
	printf("\nPSNR(R) = %f\nPSNR(G) = %f\nPSNR(B) = %f\n", PSNR(components_rgb[0], components_rgb_recovered[0]), 
	PSNR(components_rgb[1], components_rgb_recovered[1]), PSNR(components_rgb[2], components_rgb_recovered[2]));

	//8.a
	float** Cb_decimated1 = decimation1(components_ycbcr[1], 2);
	float** Cr_decimated1 = decimation1(components_ycbcr[2], 2);
	printf("\ndecimation(2) with removing of rows and columns done!\n");

	//8.b
	float** Cb_decimated2 = decimation2(components_ycbcr[1], 2);
	float** Cr_decimated2 = decimation2(components_ycbcr[2], 2);
	printf("\ndecimation(2) with average done!\n");

	//9
	std::vector<float**> components_ycbcr_recovered1;
	components_ycbcr_recovered1.push_back(components_ycbcr[0]);
	components_ycbcr_recovered1.push_back(decimation_recover(Cb_decimated1, 2));
	components_ycbcr_recovered1.push_back(decimation_recover(Cr_decimated1, 2));
	std::vector<float**> components_rgb_recovered1 = ycbcr_to_rgb(components_ycbcr_recovered1);
	printf("\nrecovery of values decimated(2) with removing of rows and columns done!\n");

	std::vector<float**> components_ycbcr_recovered2;
	components_ycbcr_recovered2.push_back(components_ycbcr[0]);
	components_ycbcr_recovered2.push_back(decimation_recover(Cb_decimated2, 2));
	components_ycbcr_recovered2.push_back(decimation_recover(Cr_decimated2, 2));
	std::vector<float**> components_rgb_recovered2 = ycbcr_to_rgb(components_ycbcr_recovered2);
	printf("\nrecovery of values decimated(2) with average done!\n");

	//10
	printf("\nPSNR[Cb1] = %f", PSNR(components_ycbcr[1], components_ycbcr_recovered1[1]));
	printf("\nPSNR[Cb2] = %f\n", PSNR(components_ycbcr[1], components_ycbcr_recovered2[1]));
	printf("\nPSNR[Cr1] = %f", PSNR(components_ycbcr[2], components_ycbcr_recovered1[2]));
	printf("\nPSNR[Cr2] = %f\n", PSNR(components_ycbcr[2], components_ycbcr_recovered2[2]));
	printf("\nPSNR[R1] = %f", PSNR(components_rgb[0], components_rgb_recovered1[0]));
	printf("\nPSNR[R2] = %f\n", PSNR(components_rgb[0], components_rgb_recovered2[0]));
	printf("\nPSNR[G1] = %f", PSNR(components_rgb[1], components_rgb_recovered1[1]));
	printf("\nPSNR[G2] = %f\n", PSNR(components_rgb[1], components_rgb_recovered2[1]));
	printf("\nPSNR[B1] = %f", PSNR(components_rgb[2], components_rgb_recovered1[2]));
	printf("\nPSNR[B2] = %f\n", PSNR(components_rgb[2], components_rgb_recovered2[2]));

	//11(8.a)
	Cb_decimated1 = decimation1(components_ycbcr[1], 4);
	Cr_decimated1 = decimation1(components_ycbcr[2], 4);
	printf("\ndecimation(4) with removing of rows and columns done!\n");

	//11(8.b)
	Cb_decimated2 = decimation2(components_ycbcr[1], 4);
	Cr_decimated2 = decimation2(components_ycbcr[2], 4);
	printf("\ndecimation(4) with average done!\n");

	//11(9)
	components_ycbcr_recovered1.clear();
	components_rgb_recovered1.clear();
	components_ycbcr_recovered1.push_back(components_ycbcr[0]);
	components_ycbcr_recovered1.push_back(decimation_recover(Cb_decimated1, 4));
	components_ycbcr_recovered1.push_back(decimation_recover(Cr_decimated1, 4));
	components_rgb_recovered1 = ycbcr_to_rgb(components_ycbcr_recovered1);
	printf("\nrecovery of values decimated(4) with removing of rows and columns done!\n");

	components_ycbcr_recovered2.clear();
	components_rgb_recovered2.clear();
	components_ycbcr_recovered2.push_back(components_ycbcr[0]);
	components_ycbcr_recovered2.push_back(decimation_recover(Cb_decimated2, 4));
	components_ycbcr_recovered2.push_back(decimation_recover(Cr_decimated2, 4));
	components_rgb_recovered2 = ycbcr_to_rgb(components_ycbcr_recovered2);
	printf("\nrecovery of values decimated(4) with average done!\n");

	//11(10)
	printf("\nPSNR[Cb1] = %f", PSNR(components_ycbcr[1], components_ycbcr_recovered1[1]));
	printf("\nPSNR[Cb2] = %f\n", PSNR(components_ycbcr[1], components_ycbcr_recovered2[1]));
	printf("\nPSNR[Cr1] = %f", PSNR(components_ycbcr[2], components_ycbcr_recovered1[2]));
	printf("\nPSNR[Cr2] = %f\n", PSNR(components_ycbcr[2], components_ycbcr_recovered2[2]));
	printf("\nPSNR[R1] = %f", PSNR(components_rgb[0], components_rgb_recovered1[0]));
	printf("\nPSNR[R2] = %f\n", PSNR(components_rgb[0], components_rgb_recovered2[0]));
	printf("\nPSNR[G1] = %f", PSNR(components_rgb[1], components_rgb_recovered1[1]));
	printf("\nPSNR[G2] = %f\n", PSNR(components_rgb[1], components_rgb_recovered2[1]));
	printf("\nPSNR[B1] = %f", PSNR(components_rgb[2], components_rgb_recovered1[2]));
	printf("\nPSNR[B2] = %f\n\n", PSNR(components_rgb[2], components_rgb_recovered2[2]));

	//12, 13
	histogram(components_rgb[0], "R", false);
	histogram(components_rgb[1], "G", false);
	histogram(components_rgb[2], "B", false);
	histogram(components_ycbcr[0], "Y", false);
	histogram(components_ycbcr[1], "Cb", false);
	histogram(components_ycbcr[2], "Cr", false);

	//14
	height_temp = height;
	width_temp = width;
	float** D_R = differential_coding(components_rgb[0]);
	float** D_G = differential_coding(components_rgb[1]);
	float** D_B = differential_coding(components_rgb[2]);
	float** D_Y = differential_coding(components_ycbcr[0]);
	float** D_Cb = differential_coding(components_ycbcr[1]);
	float** D_Cr = differential_coding(components_ycbcr[2]);

	//15, 16
	printf("\n");
	histogram(D_R, "D_R", true);
	histogram(D_G, "D_G", true);
	histogram(D_B, "D_B", true);
	histogram(D_Y, "D_Y", true);
	histogram(D_Cb, "D_Cb", true);
	histogram(D_Cr, "D_Cr", true);

	//17*(4.a)
	height = height_temp / 2;
	width = width_temp / 2;
	std::vector<float**> subframes_Y = get_subframes(components_ycbcr[0], in);

	//17*(4.b)
	height_temp = height;
	width_temp = width;
	autocorrelation(subframes_Y[0], "Y00");
	autocorrelation(subframes_Y[0], "Y01");
	autocorrelation(subframes_Y[0], "Y10");
	autocorrelation(subframes_Y[0], "Y11");
	printf("\nautocorrelation functions for Y(i,j) done!\n");

	//17*(4.c)
	printf("\ncorrelations coeffitients:\nr[Y00,Y01]=%f\n", correlation(subframes_Y[0], subframes_Y[1]));
	printf("r[Y00,Y10]=%f\n", correlation(subframes_Y[0], subframes_Y[2]));
	printf("r[Y00,Y11]=%f\n", correlation(subframes_Y[0], subframes_Y[3]));
	printf("r[Y01,Y10]=%f\n", correlation(subframes_Y[1], subframes_Y[2]));
	printf("r[Y01,Y11]=%f\n", correlation(subframes_Y[1], subframes_Y[3]));
	printf("r[Y10,Y11]=%f\n", correlation(subframes_Y[2], subframes_Y[3]));

	//17*(4.d)
	printf("\n");
	histogram(subframes_Y[0], "Y00", false);
	histogram(subframes_Y[1], "Y01", false);
	histogram(subframes_Y[2], "Y10", false);
	histogram(subframes_Y[3], "Y11", false);
	
	//17*(4.e)
	subframes_modification(subframes_Y);
	printf("\n");
	histogram(subframes_Y[0], "D_Y00", true);
	histogram(subframes_Y[1], "D_Y01", true);
	histogram(subframes_Y[2], "D_Y10", true);
	histogram(subframes_Y[3], "D_Y11", true);
	return 0;
}

void subframes_modification(std::vector<float**>& subframes)
{
	std::vector<float**> result(4);
	result[0] = subframes[0];
	for (unsigned int k = 1; k < subframes.size(); k++)
	{
		result[k] = new float*[height];
		for (int i = 0; i < height; i++)
		{
			result[k][i] = new float[width];
			for (int j = 0; j < width; j++)
			{
				result[k][i][j] = k == 1 ? subframes[k][i][j] - subframes[0][i][j] 
								: k == 2 ? subframes[k][i][j] - std::round((subframes[0][i][j] + subframes[1][i][j]) / 2)
								: subframes[k][i][j] - std::round((subframes[0][i][j] + subframes[1][i][j] + subframes[2][i][j]) / 3); 
			}
		}
	}

	subframes = result;
}

std::vector<float**> get_subframes(float** component, FILE* in)
{
	BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;

	fread(&bfh, sizeof(bfh), 1, in);
	bfh.bfSize = sizeof(bfh) + sizeof(bih) + width * height * 3;

    fread(&bih, sizeof(bih), 1, in);
    bih.biWidth = width;
    bih.biHeight = height;

    fclose(in);

    int add = 0;
    if ((bih.biWidth * 3) % 4 != 0)
    {
        add = 4 - (bih.biWidth * 3) % 4;
    }

	std::vector<float**> result;
	RGBTRIPLE rgb;

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			FILE* out;
			if (i == 0 && j == 0)
				out = fopen("Y_subframe_0_0.bmp", "wb");
			else if (i == 0 && j == 1)
				out = fopen("Y_subframe_0_1.bmp", "wb");
			else if (i == 1 && j == 0)
				out = fopen("Y_subframe_1_0.bmp", "wb");
			else
				out = fopen("Y_subframe_1_1.bmp", "wb");

			fwrite(&bfh, sizeof(bfh), 1, out);
			fwrite(&bih, sizeof(bih), 1, out);

			float** subframe = new float*[height];

			for (int y = 0; y < height; y++)
			{
				subframe[y] = new float[width];
				for (int x = 0; x < width; x++)
				{
					subframe[y][x] = component[2 * y + i][2 * x + j];
					rgb.rgbtRed = (BYTE)subframe[y][x];
					rgb.rgbtBlue = (BYTE)subframe[y][x];
					rgb.rgbtGreen = (BYTE)subframe[y][x];
					fwrite(&rgb, sizeof(rgb), 1, out);
				}
				if (add != 0)
        		{
        		    fwrite(&rgb, add, 1, out);
        		}
			}

			result.push_back(subframe);
			fclose(out);
		}
	}

	return result;
}

float** differential_coding(float** component) //r = 3
{
	height--;
	width--;
	float** result = new float*[height];

	for (int i = 0; i < height; i++)
	{
		result[i] = new float[width];
		for (int j = 0; j < width; j++)
		{
			result[i][j] = component[i + 1][j + 1] - component[i][j];
		}
	}

	return result;

}

void enthropy(std::map<float, unsigned int>& values, std::string& component_name)
{
	float enthropy = 0;
	for (std::map<float, unsigned int>::iterator i = values.begin(); i != values.end(); ++i)
	{
	 	float p = (*i).second / (float)(height * width);
	 	enthropy += p * std::log2(p);
	}
	enthropy = -enthropy;
	std::cout << "H[" << component_name << "] = " << enthropy << std::endl;
}

void histogram(float** component, std::string component_name, bool is_diff)
{
	std::map<float, unsigned int> hist;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++) 
		{
			if (hist.find(component[i][j]) == hist.end())
				hist.emplace(component[i][j], 1);
			else
			{
				unsigned int n = hist[component[i][j]] + 1;
				hist.erase(component[i][j]);
				hist.emplace(component[i][j], n);
			}
				
		}
	}

	enthropy(hist, component_name);

	if (!is_diff)
	{
		if (hist.find(0.0) == hist.end())
			hist.emplace(0.0, 0);
		if (hist.find(255.0) == hist.end())
			hist.emplace(255.0, 0);
	}

	std::ofstream hist_file;
	hist_file.open("hist" + component_name + ".csv");
	for (std::map<float, unsigned int>::iterator i = hist.begin(); i != hist.end(); ++i)
	{
		hist_file << (*i).first << "," << (*i).second << "\n";
	}
	hist_file.close();
}

float** decimation_recover(float** component, int coef)
{
	float** result = new float*[height];
	for (int i = 0; i < height; i++)
		result[i] = new float[width];

	for (int i = 0; i < height / coef; i++)
	{
		for (int j = 0; j < width / coef; j++)
		{
			for (int k = i * coef; k < i * coef + coef; k++)
				for (int l = j * coef; l < j * coef + coef; l++)
					result[k][l] = component[i][j];
		}
	}

	return result;
}

float** decimation2(float** component, int coef)
{
	float** result = new float*[height / coef];

	for (int i = 0; i < height / coef; i++)
	{
		result[i] = new float[width / coef];
		for (int j = 0; j < width / coef; j++)
		{
			result[i][j] = 0;
			for (int k = i * coef; k < i * coef + coef; k++)
				for (int l = j * coef; l < j * coef + coef; l++)
					result[i][j] += component[k][l];
			result[i][j] = std::round(result[i][j] / (coef * coef));
		}
	}

	return result;

}

float** decimation1(float** component, int coef)
{
	float** result = new float*[height / coef];
	for (int i = 0; i < height / coef; i++)
		result[i] = new float[width / coef];

	int x = 0;

	for (int i = 0; i < height; i++)
	{
		for (int j = 0, y = 0; j < width; j++)
		{
			if ((j + 1) % coef == 0 && (i + 1) % coef == 0)
			{
				result[x][y] = component[i][j];
				y++;
				if (y == width / coef)
					x++;
			}
		}
	}

	return result;
}

float PSNR(float** original, float** recovered)
{
	float psnr = 0;
	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			psnr += std::pow(original[i][j] - recovered[i][j], 2);

	psnr = 10 * (std::log10(width) + std::log10(height) + std::log10(255) + std::log10(255) - std::log10(psnr));
	return psnr;
}

std::vector<float**> ycbcr_to_rgb(std::vector<float**> ycbcr)
{
	float** R = new float*[height];
    float** G = new float*[height];
    float** B = new float*[height];

    for (int i = 0; i < height; i++)
    {
    	R[i] = new float[width];
    	G[i] = new float[width];
    	B[i] = new float[width];

    	for (int j = 0; j < width; j++)
    	{
    		R[i][j] = std::round(ycbcr[0][i][j] + 1.402 * (ycbcr[2][i][j] - 128));
    		R[i][j] = R[i][j] < 0 ? 0 : R[i][j] > 255 ? 255 : R[i][j];

    		G[i][j] = std::round(ycbcr[0][i][j] - 0.714 * (ycbcr[2][i][j] - 128) - 0.334 * (ycbcr[1][i][j] - 128));
    		G[i][j] = G[i][j] < 0 ? 0 : G[i][j] > 255 ? 255 : G[i][j];

    		B[i][j] = std::round(ycbcr[0][i][j] + 1.772 * (ycbcr[1][i][j] - 128));
    		B[i][j] = B[i][j] < 0 ? 0 : B[i][j] > 255 ? 255 : B[i][j];
    	}
    }

    std::vector<float**> result;
    result.push_back(R);
    result.push_back(G);
    result.push_back(B);

    return result;
}

std::vector<float**> rgb_to_ycbcr(std::vector<float**> rgb, FILE* in)
{
	BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;

	FILE* Y_comp = fopen("Y_comp.bmp", "wb");
	FILE* Cb_comp = fopen("Cb_comp.bmp", "wb");
	FILE* Cr_comp = fopen("Cr_comp.bmp", "wb");

	fread(&bfh, sizeof(bfh), 1, in);
    fwrite(&bfh, sizeof(bfh), 1, Y_comp);
    fwrite(&bfh, sizeof(bfh), 1, Cb_comp);
    fwrite(&bfh, sizeof(bfh), 1, Cr_comp);

    fread(&bih, sizeof(bih), 1, in);
    fwrite(&bih, sizeof(bih), 1, Y_comp);
    fwrite(&bih, sizeof(bih), 1, Cb_comp);
    fwrite(&bih, sizeof(bih), 1, Cr_comp);

    float** Y = new float*[height];
    float** Cb = new float*[height];
    float** Cr = new float*[height];

	RGBTRIPLE rgbt;

	int add = 0;
    if ((bih.biWidth * 3) % 4 != 0)
    {
        add = 4 - (bih.biWidth * 3) % 4;
    }
	for (int i = 0; i < height; i++)
	{
		Y[i] = new float[width];
		Cb[i] = new float[width];
		Cr[i] = new float[width];
		for (int j = 0; j < width; j++)
		{
			fread(&rgbt, sizeof(rgbt), 1, in);
			Y[i][j] = std::round(0.299 * rgb[0][i][j] + 0.587 * rgb[1][i][j] + 0.114 * rgb[2][i][j]);
			rgbt.rgbtRed = (BYTE)Y[i][j];
			rgbt.rgbtGreen = (BYTE)Y[i][j];
			rgbt.rgbtBlue = (BYTE)Y[i][j];
			fwrite(&rgbt, sizeof(rgbt), 1, Y_comp);

			Cb[i][j] = std::round(0.5643 * (rgb[2][i][j] - Y[i][j]) + 128);
			rgbt.rgbtRed = (BYTE)Cb[i][j];
			rgbt.rgbtGreen = (BYTE)Cb[i][j];
			rgbt.rgbtBlue = (BYTE)Cb[i][j];
			fwrite(&rgbt, sizeof(rgbt), 1, Cb_comp);

			Cr[i][j] = std::round(0.7132 * (rgb[0][i][j] - Y[i][j]) + 128);
			rgbt.rgbtRed = (BYTE)Cr[i][j];
			rgbt.rgbtGreen = (BYTE)Cr[i][j];
			rgbt.rgbtBlue = (BYTE)Cr[i][j];
			fwrite(&rgbt, sizeof(rgbt), 1, Cr_comp);
		}
		if (add != 0)
        {
            fread(&rgbt, add, 1, in);
            fwrite(&rgbt, add, 1, Y_comp);
            fwrite(&rgbt, add, 1, Cb_comp);
            fwrite(&rgbt, add, 1, Cr_comp);
        }
	}

	std::vector<float**> result;
	result.push_back(Y);
	result.push_back(Cb);
	result.push_back(Cr);

	fseek(in, 0, SEEK_SET);
	fclose(Y_comp);
	fclose(Cb_comp);
	fclose(Cr_comp);

	return result;
}

void autocorrelation(float** component, std::string component_name)
{
	std::ofstream autocorrelation_function;

	for (int y = -10, i = 1; y <= 10; y += 5, i++)
	{
		autocorrelation_function.open("autocorrelation" + component_name + "_" + std::to_string(i) + ".csv");
		for(int x = 0; x < width / 4; x++)
		{
			float** sample = get_sample(component, y, x);
			height = height - y;
			width = width - x;
			autocorrelation_function << x << "," << y << "," << correlation(component, sample) << "\n";

			for (int j = 0; j < height; j++)
				delete[] sample[j];
			delete[] sample;

			height = height_temp;
			width = width_temp;
		}
		autocorrelation_function.close();
	}	
}

float** get_sample(float** component, int y, int x)
{
	float** result = new float*[height - y];
	for (int i = 0; i < height - y; i++)
	{
		result[i] = new float[width - x];
		for (int j = 0; j < width - x; j++)
		{
			if(i + y < 0 || j + x < 0)
				result[i][j] = 0;
			else
				result[i][j] = component[i + y][j + x];
		}
	}
	return result;
}

float correlation(float** A, float** B)
{
	float M_A = 0, M_B = 0, r = 0, sigma_A = 0, sigma_B = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			if(!((height_temp < height || width_temp < width) && (i >= height_temp || j >= width_temp)))
				M_A += A[i][j];
			M_B += B[i][j];
		}
	}
	M_A /= (height * width);
	M_B /= (height * width);

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			if((height_temp < height || width_temp < width) && (i >= height_temp || j >= width_temp))
			{
				r += (0 - M_A) * (B[i][j] - M_B);
				sigma_A += std::pow(0 - M_A, 2);
			}
			else
			{
				r += (A[i][j] - M_A) * (B[i][j] - M_B);
				sigma_A += std::pow(A[i][j] - M_A, 2);
			}
			sigma_B += std::pow(B[i][j] - M_B, 2);

		}
	}
	r /= (height * width);
	sigma_A = std::sqrt(sigma_A / (height * width - 1));
	sigma_B = std::sqrt(sigma_B / (height * width - 1));
	r /= (sigma_A * sigma_B);

	return r;
}

std::vector<float**> get_rgb_components(FILE* in)
{
	BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;

    FILE* R_comp = fopen("R_comp.bmp", "wb");
    FILE* G_comp = fopen("G_comp.bmp", "wb");
    FILE* B_comp = fopen("B_comp.bmp", "wb");

    fread(&bfh, sizeof(bfh), 1, in);
    fwrite(&bfh, sizeof(bfh), 1, R_comp);
    fwrite(&bfh, sizeof(bfh), 1, G_comp);
    fwrite(&bfh, sizeof(bfh), 1, B_comp);

    fread(&bih, sizeof(bih), 1, in);
    if (bih.biBitCount != 24)
    {
        fprintf(stdout, "\nERROR! BMP24 required.\n");
        exit (1);
    }
    fwrite(&bih, sizeof(bih), 1, R_comp);
    fwrite(&bih, sizeof(bih), 1, G_comp);
    fwrite(&bih, sizeof(bih), 1, B_comp);


    height = bih.biHeight;
    width = bih.biWidth;
    height_temp = height;
    width_temp = width;

    float** R = new float*[height];
    float** G = new float*[height];
    float** B = new float*[height];

    RGBTRIPLE rgb;

    int add = 0;
    if ((bih.biWidth * 3) % 4 != 0) 
        add = 4 - (bih.biWidth * 3) % 4;

    for (int i = 0; i < bih.biHeight; i++)
    {
    	R[i] = new float[bih.biWidth];
    	G[i] = new float[bih.biWidth];
    	B[i] = new float[bih.biWidth];

        for (int j = 0; j < bih.biWidth; j++)
        {
            fread(&rgb, sizeof(rgb), 1, in);

            RGBTRIPLE rgb_temp = rgb;
            rgb_temp.rgbtGreen = 0;
            rgb_temp.rgbtBlue = 0;
            R[i][j] = rgb_temp.rgbtRed;
            fwrite(&rgb_temp, sizeof(rgb_temp), 1, R_comp);

            rgb_temp = rgb;
            rgb_temp.rgbtRed = 0;
            rgb_temp.rgbtBlue = 0;
            G[i][j] = rgb_temp.rgbtGreen;
            fwrite(&rgb_temp, sizeof(rgb_temp), 1, G_comp);

            rgb_temp = rgb;
            rgb_temp.rgbtRed = 0;
            rgb_temp.rgbtGreen = 0;
            B[i][j] = rgb_temp.rgbtBlue;
            fwrite(&rgb_temp, sizeof(rgb_temp), 1, B_comp);
        }
        if (add != 0)
        {
            fread(&rgb, add, 1, in);
            fwrite(&rgb, add, 1, R_comp);
            fwrite(&rgb, add, 1, G_comp);
            fwrite(&rgb, add, 1, B_comp);
        }
    }

    fseek(in, 0, SEEK_SET);
    fclose(R_comp);
    fclose(G_comp);
    fclose(B_comp);

    std::vector<float**> result;
    result.push_back(R);
    result.push_back(G);
    result.push_back(B);

    return result;
}