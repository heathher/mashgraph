#include "align.h"
#include <string>

#include <climits>
#include <cmath>
using std::string;
using std::cout;
using std::endl;

//find abs(n)
int abs(int n){
    if (n < 0)
        return -n;
    else 
        return n;
}

//find max(a, b)
int max(int a, int b){
    if (a > b)
        return a;
    else 
        return b;
}

//find min(a, b)
int min(int a, int b){
    if (a < b)
        return a;
    else 
        return b;
}

//find MSE
unsigned long long metrika(const Image &im1, const Image &im2){
    uint r1 = 0, r2 = 0;
    unsigned long long sum = 0;
    for (uint i = 0; i < im1.n_rows; i++){
        for (uint j = 0; j < im1.n_cols; j++){
            r1 = std::get<0>(im1(i,j));
            r2 = std::get<0>(im2(i,j));
            sum += pow(r1 < r2 ? r2 - r1 : r1 - r2 , 2);
        }
    } 
    return (sum/(im1.n_rows*im1.n_cols));
}

//get result image 
//
Image getresImage(const Image &redImage, const Image &greenImage, const Image &blueImage, std::pair<int, int> offset1, 
                    std::pair<int, int> offset2, int shift){

    int row = redImage.n_rows;
    int col = redImage.n_cols;
    uint r, g, b;
    Image resImage(row + 2*shift, col + 2*shift);
    for (int i = shift; i < row - shift; ++i){
        for (int j = shift; j < col - shift ; ++j){
            r = std::get<0>(redImage(i, j));
            g = std::get<1>(greenImage(i + offset1.first, j + offset1.second));
            b = std::get<2>(blueImage(i + offset2.first, j + offset2.second));
            resImage(i, j) = std::make_tuple(r, g, b);
        }
    }
    resImage = resImage.submatrix(shift, shift, row - 2*shift, col - 2*shift);
    return resImage;
}

//find offset of image 
// piramida if row > 500 || col > 500
std::pair<int, int> getOffset(const Image &im1, const Image &im2, int shift){
    std::pair<int, int> offset;
    int row = im1.n_rows;
    int col = im1.n_cols;
    unsigned long long minsum = ULLONG_MAX;
    // if (row > 500 || col > 500){
    //     // halve 
    //     Image n_im1 = resize(im1, 0.5);
    //     Image n_im2 = resize(im2, 0.5);
    //     //get offset after change image
    //     //recursion
    //     offset = getOffset(n_im1, n_im2, shift, row/2, col/2);
    //     offset.first *= 2;
    //     offset.second *= 2; 
    //     for (int k = offset.first - 3; k < offset.first + 3; k++){
    //         for (int l = offset.second - 3; l < offset.second + 3; l++){
    //             Image tmp_image1 = im1.submatrix(min(0, -k), min(0, -l), row - abs(k) - 1, col - abs(l) - 1);
    //             Image tmp_image2 = im2.submatrix(max(0, k), max(0, l), row - abs(k) - 1, col - abs(l) - 1);
    //             unsigned long long sum = metrika(tmp_image1, tmp_image2, shift);
    //             if (sum < minsum){
    //                 minsum = sum;
    //                 offset.first = k;
    //                 offset.second = l;
    //             }
    //         }
    //     }
    // }
    // else 
    // {
        //if col < 500 and row < 500
        for (int k = -shift; k < shift; ++k){
            for (int l = -shift; l < shift; ++l){
                Image tmp_image1 = im1.submatrix(row*0.05, col*0.05, row*0.9, col*0.9);
                Image tmp_image2 = im2.submatrix(row*0.05+k, col*0.05+l, row*0.9, col*0.9);
                double sum = metrika(tmp_image1, tmp_image2);
                if (sum < minsum){
                    minsum = sum;
                    offset.first = k;
                    offset.second = l;
                }
            }
        }
    //}
    
    return offset;
}



Image align(const Image &src, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, 
            bool isInterp, bool isSubpixel, double subScale)
{
    int shift = 15;
    Image srcImage = src;
    uint row = srcImage.n_rows/3, col = srcImage.n_cols;

    //--subpixel increase
    if (isSubpixel){
        srcImage = subpixelBig(src, subScale);
        row = srcImage.n_rows/3;
        col = srcImage.n_cols;
    }

    // cut to red, green, blue
    Image 
        blueImage = srcImage.submatrix(0, 0, row, col),
        greenImage = srcImage.submatrix(row, 0, row, col),
        redImage = srcImage.submatrix(2*row, 0, row, col);

    //find offset between red and green -> piramida
    std::pair<int, int> offset1 = getOffset(redImage, greenImage, shift);

    //find offset between red and blue -> piramida
    std::pair<int, int> offset2 = getOffset(redImage, blueImage, shift);

    //get result image 
    Image resImage = getresImage(redImage, greenImage, blueImage, offset1, offset2, shift);

    //--subpixel reduce
    if (isSubpixel){
        resImage = subpixelSmall(resImage, subScale);
    }

    //postprocess function
    if (isPostprocessing){
        if (postprocessingType == "--gray-world")
            resImage = gray_world(resImage);
        else if (postprocessingType == "--unsharp"){
            if (isMirror){
                resImage = mirror(resImage, 1);
                resImage = unsharp(resImage);
                resImage = disMirror(resImage, 1);
            }
            else resImage = unsharp(resImage);
        }
        else if (postprocessingType == "--autocontrast")
            resImage = autocontrast(resImage, fraction);
    }
    return resImage;
}

//check a in [0, 255]
int border_check(int a){
    if (a < 0)
        return 0;
    else if (a > 255)
        return 255;
    else return a;
}

Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};
    return custom(src_image, kernel);
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    return custom(src_image, kernel);
}

//--subpixel function
Image subpixelSmall(Image src_image, double subScale){
    int row = src_image.n_rows / subScale, col = src_image.n_cols / subScale; 
    Image res_image(row, col);
    for (int i = 0; i < row; ++i){
        for (int j = 0; j < col; ++j){
            res_image(i, j) = src_image(i * subScale, j * subScale); 
        }
    }
    return res_image;
}

//--subpixel function
Image subpixelBig(Image src_image, double subScale){
    int row = src_image.n_rows * subScale, col = src_image.n_cols * subScale;
    Image res_image(row, col);
    for (int i = 0; i < row; ++i){
        for (int j = 0; j < col; ++j){
            res_image(i, j) = src_image(i / subScale,j / subScale);
        }
    }
    return res_image;
}

//--unsharp
Image unsharp(Image src_image) {
    Matrix<double> kernel = {{-1.0/6, -2.0/3, -1.0/6},
                             {-2.0/3, 13.0/3, -2.0/3},
                             {-1.0/6, -2.0/3, -1.0/6}};
    double r, g, b;
    Image tmp_image = src_image.deep_copy(); 
    for (uint i = 1; i < src_image.n_rows - 1; i++){
        for (uint j = 1; j < src_image.n_cols - 1; j++){
            double sum_r = 0, sum_g = 0, sum_b = 0;
            for (int k = -1; k < 2; k++){
                for (int l = -1; l < 2; l++){
                    std::tie(r, g, b) = tmp_image(i + k, j + l);
                    sum_r += r * kernel(k + 1, l + 1);
                    sum_g += g * kernel(k + 1, l + 1);
                    sum_b += b * kernel(k + 1, l + 1);
                }
            }
            src_image(i, j) = std::make_tuple(border_check(sum_r), border_check(sum_g), border_check(sum_b));
        } 
    }
    return src_image;
}

//cut border after mirror()
Image disMirror(Image src_image, uint offset_size){
    return src_image.submatrix(offset_size, offset_size, src_image.n_rows - 2*offset_size, src_image.n_cols - 2*offset_size);
}

//--mirror
Image mirror(Image src_image, uint offset_size){
    uint row = src_image.n_rows, col = src_image.n_cols;
    Image res_image(row + 2 * offset_size + 1, col + 2 * offset_size + 1 );
    
    for (uint i = 0; i < row; ++i){
        for (uint j = 0; j < col; ++j){
            res_image(i + offset_size, j + offset_size) = src_image(i, j);
        }
    }
    
    //left and right border
    for (uint i = 0; i < row; ++i){
        for (uint j = 0; j < offset_size; ++j){
            res_image(i + offset_size, offset_size - j - 1) = src_image(i, j);
            res_image(i + offset_size, offset_size + col + j - 1) = src_image(i, col - j - 1);
        }
    }

    //upper and lower border
    for (uint i = 0; i < offset_size; ++i){
        for (uint j = 0; j < col + 2 * offset_size; ++j){
            res_image(offset_size - i - 1, j) = res_image(i + offset_size, j);
            res_image(offset_size + row + i - 1, j) = res_image(row - i - 1 + offset_size, j);
        }
    }
    return res_image;
}

Image gray_world(Image src_image) {
    unsigned long long sum_r = 0, sum_b = 0, sum_g = 0, s = 0;
    uint new_r, new_g, new_b, r, g, b;
    uint row = src_image.n_rows, col = src_image.n_cols;

    //find sum r, g, b pixels
    for (uint i = 0; i < row; ++i)
        for (uint j = 0; j < col; j++){
            std::tie(r, g, b) = src_image(i, j);
            sum_r += r;
            sum_g += g;
            sum_b += b;
        }
    sum_r /= row * col;
    sum_g /= row * col;
    sum_b /= row * col;

    //find S
    s = (sum_r + sum_b + sum_g)/3;

    //make gray picture
    for (uint i = 0; i < row; i++){
        for (uint j = 0; j < col; j++){
            std::tie(r, g, b) = src_image(i, j);
            new_r = r * s / sum_r;
            new_g = g * s / sum_g;
            new_b = b * s / sum_b;
            src_image(i, j) = std::make_tuple(border_check(new_r),border_check(new_g), border_check(new_b));
        }
    }
    return src_image;
}


Image resize(Image src_image, double scale) {
    uint row = (src_image.n_rows - 1) * scale , col = (src_image.n_cols - 1) * scale ;
    Image res_image(row, col);
    uint r1, g1, b1, r2, g2, b2;

    for (uint i = 0; i < row; ++i){
        for (uint j = 0; j < col; ++j){
            //билинейное интерполирование
            double new_i = i / scale;
            double new_j = j / scale;
            int y = (j/scale);
            int x = (i/scale);
            double alpha = new_i - x;
            std::tie(r1, g1, b1) = src_image(x + 1, y);
            std::tie(r2, g2, b2) = src_image(x, y);
            std::tuple<double, double, double> value1 = std::make_tuple(r1 * alpha + r2 * (1 - alpha), 
                                                                        g1 * alpha + g2 * (1 - alpha), 
                                                                        b1 * alpha + b2 * (1 - alpha));
            std::tie(r1, g1, b1) = src_image(x + 1, y + 1);
            std::tie(r2, g2, b2) = src_image(x, y + 1);
            std::tuple<double, double, double> value2 = std::make_tuple(r1 * alpha + r2 * (1 - alpha), 
                                                                        g1 * alpha + g2 * (1 - alpha), 
                                                                        b1 * alpha + b2 * (1 - alpha));
            double beta = new_j - y;
            std::tie(r1, g1, b1) = value2;
            std::tie(r2, g2, b2) = value1;
            res_image(i, j) = std::make_tuple(r1 * beta + r2 * (1 - beta), 
                                                g1 * beta + g2 * (1 - beta),
                                                b1 * beta + b2 * (1 - beta)); 
        }
    }
    return res_image;
}

Image custom(Image src_image, Matrix<double> kernel) {
    return src_image;
}

Image autocontrast(Image src_image, double fraction) {
    std::vector<double> y;
    uint r, g, b, new_r, new_g, new_b;
    double tmp_y;
    uint row = src_image.n_rows, col = src_image.n_cols;
    uint fraction_step = fraction * row * col;

    //make array of brightness - y
    for (uint i = 0; i < row; ++i){
        for (uint j = 0; j < col; ++j){
            std::tie(r, g, b) = src_image(i, j);
            tmp_y = r * 0.2125 + g * 0.7154 + b * 0.0721;
            y.push_back(tmp_y);
        }
    }

    //sort y
    std::sort(y.begin(), y.end());

    //find min and max of y 
    double y_min = y[fraction_step], y_max = y[row * col - 1 - fraction_step];
    for (uint i = 0; i < row; ++i){
        for (uint j = 0; j < col; ++j){
            std::tie(r, g, b) = src_image(i, j);
            new_r = (r - y_min) * 255.0 / (y_max - y[0]);
            new_g = (g - y_min) * 255.0 / (y_max - y[0]);
            new_b = (b - y_min) * 255.0 / (y_max - y[0]); 
            src_image(i, j) = std::make_tuple(border_check(new_r), border_check(new_g), border_check(new_b));
        }    
    }
    return src_image;
}

Image gaussian(Image src_image, double sigma, int radius)  {
    return src_image;
}

Image gaussian_separable(Image src_image, double sigma, int radius) {
    return src_image;
}

Image median(Image src_image, int radius) {                                       
    uint size = 2 * radius + 1;                                                       
    uint row = src_image.n_rows, col = src_image.n_cols;
    uint r = 0, g = 0, b = 0, r_med, g_med, b_med;
    Image tmp_image = src_image.deep_copy();
    for (uint i = 0; i < row - size; ++i){
        for (uint j = 0; j < col - size; ++j){
            std::vector<int> r_hist(256, 0), g_hist(256, 0), b_hist(256, 0); 
            // take histogram
            for (uint k = 0; k < size; ++k){                                           
                for (uint l = 0; l < size; ++l){
                    std::tie(r, g, b) = tmp_image(i + k, j + l);
                    r_hist[r]++;
                    g_hist[g]++;
                    b_hist[b]++;
                }
            }
            uint r_pixels_count = 0;
            uint g_pixels_count = 0;
            uint b_pixels_count = 0;
            uint check_size = size * size / 2;
            //find r median
            for (uint k = 0; k < 256; ++k){
                r_pixels_count += r_hist[k];
                if (r_pixels_count >= check_size){
                    r_med = k;
                    break;
                }
            }
            //find g median
            for (uint k = 0; k < 256; ++k){
                g_pixels_count += g_hist[k];
                if (g_pixels_count >= check_size){
                    g_med = k;
                    break;
                }
            }
            //find b median
            for (uint k = 0; k < 256; ++k){
                b_pixels_count += b_hist[k];
                if (b_pixels_count >= check_size){
                    b_med = k;
                    break;
                }
            }
            // put result in image
            src_image(i + radius, j + radius) = std::make_tuple(r_med, g_med, b_med);
        }
    }
    return src_image;
}

Image median_linear(Image src_image, int radius) {
    uint size = 2 * radius + 1;                                            
    uint row = src_image.n_rows, col = src_image.n_cols;          
    uint r, g, b, r_med, g_med, b_med;
    Image tmp_image = src_image.deep_copy();
    for (uint i = 0; i < row - size; ++i){
        std::vector<uint> r_hist(256, 0), g_hist(256, 0), b_hist(256, 0);
        // get first histogram in a row
        for (uint k = 0; k < size; ++k){                      
            for (uint l = 0; l < size; ++l){
                std::tie(r, g, b) = tmp_image(k + i, l);
                r_hist[r]++;
                g_hist[g]++;
                b_hist[b]++;
            }
        }
        // looking in a col
        for (uint j = 0; j < col - size; ++j){
            uint r_pixels_count = 0;
            uint g_pixels_count = 0;
            uint b_pixels_count = 0;
            uint check_size = size * size / 2;
            // find r median
            for (uint k = 0; k < 256; ++k){
                r_pixels_count += r_hist[k];
                if (r_pixels_count >= check_size){
                    r_med = k;
                    break;
                }
            }
            //find g median
            for (uint k = 0; k < 256; ++k){
                g_pixels_count += g_hist[k];
                if (g_pixels_count >= check_size){
                    g_med = k;
                    break;
                }
            }
            //find b median
            for (uint k = 0; k < 256; ++k){
                b_pixels_count += b_hist[k];
                if (b_pixels_count >= check_size){
                    b_med = k;
                    break;
                }
            }
            //put result in image
            src_image(i + radius, j + radius) = std::make_tuple(r_med, g_med, b_med);
            // delete unusable row from histogram
            for (uint k = 0; k < size; ++k){                                            
                std::tie(r, g, b) = tmp_image(i + k,j);
                r_hist[r]--;
                g_hist[g]--;
                b_hist[b]--;
            }
            //add new row to histogram
            for (uint k = 0; k < size; ++k){
                std::tie(r, g, b) = tmp_image(i + k, j + size);
                r_hist[r]++;
                g_hist[g]++;
                b_hist[b]++;
            }
        }

    }
    return src_image;
}

Image median_const(Image src_image, int radius) {
    return src_image;
}

Image canny(Image src_image, int threshold1, int threshold2) {
    return src_image;
}
