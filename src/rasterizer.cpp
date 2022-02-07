#include "rasterizer.h"

using namespace std;

namespace CGL {

RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
                             size_t width, size_t height,
                             unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;
    
    sample_buffer.resize(width * height * sample_rate, Color::White);
}

// Used by rasterize_point and rasterize_line
void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
    for (size_t dy=0; dy<supersampling_size; dy++) {
        for (size_t dx=0; dx<supersampling_size; dx++) {
            sample_buffer[(y * width + x) * sample_rate + dy * supersampling_size + dx] = c;
        }
    }
}

// Rasterize a point: simple example to help you start familiarizing
// yourself with the starter code.
//
void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);
    
    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;
    
    fill_pixel(sx, sy, color);
    return;
}

// Rasterize a line.
void RasterizerImp::rasterize_line(float x0, float y0,
                                   float x1, float y1,
                                   Color color) {
    if (x0 > x1) {
        swap(x0, x1); swap(y0, y1);
    }
    
    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
        dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
        dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }
    
    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
        rasterize_point(pt[0], pt[1], color);
        pt[0] += dpt[0]; pt[1] += dpt[1];
    }
}

// Rasterize a triangle.
void RasterizerImp::rasterize_triangle(float x0, float y0,
                                       float x1, float y1,
                                       float x2, float y2,
                                       Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    // Normal vector of three lines: n = (-(y1-y0), (x1-x0))
    Vector2D n0(-(y1-y0), x1-x0), n1(-(y2-y1), x2-x1), n2(-(y0-y2), x0-x2);
    
    // Bounding box of the triangle
    size_t x_min = min(min(floor(x0), floor(x1)), min(floor(x2), floor((float) width)));
    size_t x_max = max(max(ceil(x0), ceil(x1)), max(ceil(x2), ceil((float) 0)));
    size_t y_min = min(min(floor(y0), floor(y1)), min(floor(y2), floor((float) height)));
    size_t y_max = max(max(ceil(y0), ceil(y1)), max(ceil(y2), ceil((float) 0)));
    
    // Loop through x and y in the bounding box
    float step_size = 1.0 / supersampling_size;
    
    for (size_t x=x_min; x<x_max; x++) {
        for (size_t y=y_min; y<y_max; y++) {
            
            for (size_t dx=0; dx<supersampling_size; dx++) {
                for (size_t dy=0; dy<supersampling_size; dy++) {
                    float x_center = x + step_size * 0.5 + dx * step_size;
                    float y_center = y + step_size * 0.5 + dy * step_size;
                    Vector2D v0(x_center-x0, y_center-y0), v1(x_center-x1, y_center-y1), v2(x_center-x2, y_center-y2);
                    // Test whether sample point is in our triangle
                    float v0_n0 = dot(v0, n0), v1_n1 = dot(v1, n1), v2_n2 = dot(v2, n2);
                    if ((v0_n0 >= 0 && v1_n1 >= 0 && v2_n2 >= 0) ||
                        (v0_n0 <= 0 && v1_n1 <= 0 && v2_n2 <= 0)) {
                        sample_buffer[(y*width+x)*sample_rate + dy*supersampling_size+dx] = color;
                    }
                }
            }
        }
    }
    // TODO: Task 1: Extra credit: faster triangle rasterization
    // TODO: Task 2: Update to implement super-sampled rasterization
    // If you choose not to use the full width * height * sample_rate buffer space,
    // you will have to maintain another data structure to keep track of which pixel
    // has been supersampled how many times. Otherwise you will get a gray line between
    // two adjacent black triangles when sample rate is high.
    
    // I choose to use the full width * height * sample_rate buffer space because it
    // saves a bunch of if-elses, and I only need to perform one average operation
    // at the very end
}


void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
                                                          float x1, float y1, Color c1,
                                                          float x2, float y2, Color c2)
{
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
    // Normal vector of three lines: n = (-(y1-y0), (x1-x0))
    Vector2D n0(-(y1-y0), x1-x0), n1(-(y2-y1), x2-x1), n2(-(y0-y2), x0-x2);
    // Bounding box of the triangle
    size_t x_min = min(min(floor(x0), floor(x1)), min(floor(x2), floor((float) width)));
    size_t x_max = max(max(ceil(x0), ceil(x1)), max(ceil(x2), ceil((float) 0)));
    size_t y_min = min(min(floor(y0), floor(y1)), min(floor(y2), floor((float) height)));
    size_t y_max = max(max(ceil(y0), ceil(y1)), max(ceil(y2), ceil((float) 0)));
    // Ranges of alpha and bete
    float alpha_max = dot(Vector2D(x0 - x1, y0 - y1), n1);
    float beta_max = dot(Vector2D(x1 - x2, y1 - y2), n2);
    // Loop through x and y in the bounding box
    for (size_t x=x_min; x<x_max; x++) {
        for (size_t y=y_min; y<y_max; y++) {
            float x_center = x + 0.5, y_center = y + 0.5;
            Vector2D v0(x_center-x0, y_center-y0), v1(x_center-x1, y_center-y1), v2(x_center-x2, y_center-y2);
            // Test whether sample point is in our triangle
            float v0_n0 = dot(v0, n0), v1_n1 = dot(v1, n1), v2_n2 = dot(v2, n2);
            if ((v0_n0 >= 0 && v1_n1 >= 0 && v2_n2 >= 0) ||
                (v0_n0 <= 0 && v1_n1 <= 0 && v2_n2 <= 0)) {
                float alpha = v1_n1 / alpha_max;
                float beta = v2_n2 / beta_max;
                float gamma = 1 - alpha - beta;
                fill_pixel(x, y, alpha * c0 + beta * c1 + gamma * c2);
            }
        }
    }
}


void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
                                                float x1, float y1, float u1, float v1,
                                                float x2, float y2, float u2, float v2,
                                                Texture& tex)
{
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
    // uv vectors of point 0, 1, 2
    Vector2D uv0(u0, v0), uv1(u1, v1), uv2(u2, v2);
    // Normal vector of three lines: n = (-(y1-y0), (x1-x0))
    Vector2D n0(-(y1-y0), x1-x0), n1(-(y2-y1), x2-x1), n2(-(y0-y2), x0-x2);
    
    // Bounding box of the triangle
    size_t x_min = min(min(floor(x0), floor(x1)), min(floor(x2), floor((float) width)));
    size_t x_max = max(max(ceil(x0), ceil(x1)), max(ceil(x2), ceil((float) 0)));
    size_t y_min = min(min(floor(y0), floor(y1)), min(floor(y2), floor((float) height)));
    size_t y_max = max(max(ceil(y0), ceil(y1)), max(ceil(y2), ceil((float) 0)));
    // Ranges of alpha and bete
    float alpha_max = dot(Vector2D(x0 - x1, y0 - y1), n1);
    float beta_max = dot(Vector2D(x1 - x2, y1 - y2), n2);
    // Loop through x and y in the bounding box
    float step_size = 1.0 / supersampling_size;
    
    for (size_t x=x_min; x<x_max; x++) {
        for (size_t y=y_min; y<y_max; y++) {
            float pixel_x = x + 0.5, pixel_y = y + 0.5;
            // Setting up params
            SampleParams params;
            
            Vector2D pixel_v1(pixel_x-x1, pixel_y-y1), pixel_v2(pixel_x-x2, pixel_y-y2);
            float alpha = dot(pixel_v1, n1) / alpha_max;
            float beta = dot(pixel_v2, n2) / beta_max;
            float gamma = 1 - alpha - beta;
            params.p_uv = alpha * uv0 + beta * uv1 + gamma * uv2;
            
            // barycentric uv for (x, y+1)
            alpha = dot(pixel_v1 + Vector2D(0, 1), n1) / alpha_max;
            beta = dot(pixel_v2 + Vector2D(0, 1), n2) / beta_max;
            gamma = 1 - alpha - beta;
            params.p_dy_uv = alpha * uv0 + beta * uv1 + gamma * uv2;
            
            // barycentric uv for (x+1, y)
            alpha = dot(pixel_v1 + Vector2D(1, 0), n1) / alpha_max;
            beta = dot(pixel_v2 + Vector2D(1, 0), n2) / beta_max;
            gamma = 1 - alpha - beta;
            params.p_dx_uv = alpha * uv0 + beta * uv1 + gamma * uv2;
            
            params.psm = this->psm;
            params.lsm = this->lsm;
            // Sample the texture
            Color color = tex.sample(params);
            
            for (size_t dx=0; dx<supersampling_size; dx++) {
                for (size_t dy=0; dy<supersampling_size; dy++) {
                    float x_center = x + step_size * 0.5 + dx * step_size;
                    float y_center = y + step_size * 0.5 + dy * step_size;
                    Vector2D v0(x_center-x0, y_center-y0), v1(x_center-x1, y_center-y1), v2(x_center-x2, y_center-y2);
                    // Test whether sample point is in our triangle
                    float v0_n0 = dot(v0, n0), v1_n1 = dot(v1, n1), v2_n2 = dot(v2, n2);
                    if ((v0_n0 >= 0 && v1_n1 >= 0 && v2_n2 >= 0) ||
                        (v0_n0 <= 0 && v1_n1 <= 0 && v2_n2 <= 0)) {
                        sample_buffer[(y*width+x)*sample_rate + dy*supersampling_size+dx] = color;
                    }
                }
            }
        }
    }
}


void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support
    
    this->sample_rate = rate;
    this->supersampling_size = round(sqrt(rate));
    
    this->sample_buffer.resize(width * height * sample_rate, Color::White);
}


void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
                                           size_t width, size_t height)
{
    // TODO: Task 2: You may want to update this function for supersampling support
    
    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;
    
    
    this->sample_buffer.resize(width * height, Color::White);
}


void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
}


// This function is called at the end of rasterizing all elements of the
// SVG file.  If you use a supersample buffer to rasterize SVG elements
// for antialising, you could use this call to fill the target framebuffer
// pixels from the supersample buffer data.
//
void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support
    
    size_t supersampling_size = round(sqrt(this->sample_rate));
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            Color col;
            for (size_t dy=0; dy<supersampling_size; dy++) {
                for (size_t dx=0; dx<supersampling_size; dx++) {
                    col += sample_buffer[(y * width + x)*sample_rate + dy*supersampling_size+dx];
                }
            }
            col *= 1.0 / sample_rate;
            
            for (int k = 0; k < 3; ++k) {
                this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
            }
        }
    }
    
}

Rasterizer::~Rasterizer() { }


}// CGL
