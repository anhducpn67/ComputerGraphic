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
        for (int k = 0; k < sample_rate; k++) {
            sample_buffer[(y * width + x) * sample_rate + k] = c;
        }
    }

    // Rasterize a point: simple example to help you start familiarizing
    // yourself with the starter code.
    //
    void RasterizerImp::rasterize_point(float x, float y, Color color) {
        // fill in the nearest pixel
        int sx = (int) floor(x);
        int sy = (int) floor(y);
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
            swap(x0, x1);
            swap(y0, y1);
        }
        float pt[] = {x0, y0};
        float m = (y1 - y0) / (x1 - x0);
        float dpt[] = {1, m};
        int steep = abs(m) > 1;
        if (steep) {
            dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
            dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
        }

        while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
            rasterize_point(pt[0], pt[1], color);
            pt[0] += dpt[0];
            pt[1] += dpt[1];
        }
    }

    bool inside_triangle(float x, float y,
                         float x0, float y0,
                         float x1, float y1,
                         float x2, float y2) {
        double L1 = -(x - x0) * (y1 - y0) + (y - y0) * (x1 - x0);
        double L2 = -(x - x1) * (y2 - y1) + (y - y1) * (x2 - x1);
        double L3 = -(x - x2) * (y0 - y2) + (y - y2) * (x0 - x2);
        return ((L1 >= 0.0f) && (L2 >= 0.0f) && (L3 >= 0.0f));
    }

    typedef std::pair<float, float> ii;

    bool ccw(ii O, ii A, ii B) {
        A.first -= O.first;
        A.second -= O.second;
        B.first -= O.first;
        B.second -= O.second;
        return A.first * B.second > A.second * B.first;
    }

    // Rasterize a triangle.
    void RasterizerImp::rasterize_triangle(float x0, float y0,
                                           float x1, float y1,
                                           float x2, float y2,
                                           Color color) {
        if (!ccw(ii(x0, y0), ii(x1, y1), ii(x2, y2))) {
            swap(x1, x2);
            swap(y1, y2);
        }
        int min_x = (int) min(min(x0, x1), x2) - 1;
        int max_x = (int) max(max(x0, x1), x2) + 1;
        int min_y = (int) min(min(y0, y1), y2) - 1;
        int max_y = (int) max(max(y0, y1), y2) + 1;
        // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
        //    for (int x = max(0, min_x); x < min((int)width, max_x); x++) {
        //        for (int y = max(0, min_y); y < min((int)height, max_y); y++) {
        //            float center_x = 1.0*x + 0.5;
        //            float center_y = 1.0*y + 0.5;
        //            if (inside_triangle(center_x, center_y, x0, y0, x1, y1, x2, y2)) {
        //                fill_pixel(x, y, color);
        //            }
        //        }
        //    }
        // TODO: Task 2: Update to implement super-sampled rasterization
        int scale = sqrt(sample_rate);
        for (int x = max(0, min_x); x < min((int) width, max_x); x++) {
            for (int y = max(0, min_y); y < min((int) height, max_y); y++) {
                for (int i_x = 0; i_x < scale; i_x++)
                    for (int i_y = 0; i_y < scale; i_y++) {
                        float center_x = x + i_x * (1.0 / scale) + 1.0 / (2 * scale);
                        float center_y = y + i_y * (1.0 / scale) + 1.0 / (2 * scale);
                        if (inside_triangle(center_x, center_y, x0, y0, x1, y1, x2, y2)) {
                            int k = i_x * scale + i_y;
                            sample_buffer[(y * width + x) * sample_rate + k] = color;
                        }
                    }
            }
        }
    }

    float L(Vector2D origin, Vector2D P, Vector2D Q) {
        return -(origin.x - P.x) * (Q.y - P.y) + (origin.y - P.y) * (Q.x - P.x);
    }

    void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
                                                              float x1, float y1, Color c1,
                                                              float x2, float y2, Color c2) {
        // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
        // Hint: You can reuse code from rasterize_triangle
        if (!ccw(ii(x0, y0), ii(x1, y1), ii(x2, y2))) {
            swap(x1, x2);
            swap(y1, y2);
        }
        int min_x = (int) min(min(x0, x1), x2) - 1;
        int max_x = (int) max(max(x0, x1), x2) + 1;
        int min_y = (int) min(min(y0, y1), y2) - 1;
        int max_y = (int) max(max(y0, y1), y2) + 1;
        Vector2D A(x0, y0), B(x1, y1), C(x2, y2);
        int scale = sqrt(sample_rate);
        for (int x = max(0, min_x); x < min((int) width, max_x); x++) {
            for (int y = max(0, min_y); y < min((int) height, max_y); y++) {
                for (int i_x = 0; i_x < scale; i_x++)
                    for (int i_y = 0; i_y < scale; i_y++) {
                        float center_x = x + i_x * (1.0 / scale) + 1.0 / (2 * scale);
                        float center_y = y + i_y * (1.0 / scale) + 1.0 / (2 * scale);
                        Vector2D O(center_x, center_y);
                        if (inside_triangle(center_x, center_y, x0, y0, x1, y1, x2, y2)) {
                            int k = i_x * scale + i_y;
                            float alpha = L(O, B, C) / L(A, B, C);
                            float beta = L(O, A, C) / L(B, A, C);
                            float gamma = 1 - alpha - beta;
                            sample_buffer[(y * width + x) * sample_rate + k] = alpha * c0 + beta * c1 + gamma * c2;
                        }
                    }
            }
        }
    }

    void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
                                                    float x1, float y1, float u1, float v1,
                                                    float x2, float y2, float u2, float v2,
                                                    Texture &tex) {
//        // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
        if (!ccw(ii(x0, y0), ii(x1, y1), ii(x2, y2))) {
            swap(x1, x2);
            swap(u1, u2);
            swap(y1, y2);
            swap(v1, v2);
        }
        int min_x = (int) min(min(x0, x1), x2) - 1;
        int max_x = (int) max(max(x0, x1), x2) + 1;
        int min_y = (int) min(min(y0, y1), y2) - 1;
        int max_y = (int) max(max(y0, y1), y2) + 1;
        Vector2D A(x0, y0), B(x1, y1), C(x2, y2);
        Vector2D A_texture(u0, v0), B_texture(u1, v1), C_texture(u2, v2);
        int scale = sqrt(sample_rate);
        for (int x = max(0, min_x); x < min((int) width, max_x); x++) {
            for (int y = max(0, min_y); y < min((int) height, max_y); y++) {
                for (int i_x = 0; i_x < scale; i_x++)
                    for (int i_y = 0; i_y < scale; i_y++) {
                        float center_x = x + i_x * (1.0 / scale) + 1.0 / (2 * scale);
                        float center_y = y + i_y * (1.0 / scale) + 1.0 / (2 * scale);
                        Vector2D O(center_x, center_y);
                        if (inside_triangle(center_x, center_y, x0, y0, x1, y1, x2, y2)) {
                            int k = i_x * scale + i_y;
                            float alpha = L(O, B, C) / L(A, B, C);
                            float beta = L(O, A, C) / L(B, A, C);
                            float gamma = 1 - alpha - beta;
                            Vector2D O_texture = alpha * A_texture + beta * B_texture + gamma * C_texture;
                            sample_buffer[(y * width + x) * sample_rate + k] = tex.sample_nearest(O_texture);
                        }
                    }
            }
        }
        // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
        // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle

    }

    void RasterizerImp::set_sample_rate(unsigned int rate) {
        // TODO: Task 2: You may want to update this function for supersampling support (OK)
        this->sample_rate = rate;
        this->sample_buffer.resize(width * height * sample_rate, Color::White);
    }

    void RasterizerImp::set_framebuffer_target(unsigned char *rgb_framebuffer,
                                               size_t width, size_t height) {
        // TODO: Task 2: You may want to update this function for supersampling support (OK)
        this->width = width;
        this->height = height;
        this->rgb_framebuffer_target = rgb_framebuffer;
        this->sample_buffer.resize(width * height * sample_rate, Color::White);
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
        // Task 1
        //      for (int x = 0; x < width; ++x) {
        //          for (int y = 0; y < height; ++y) {
        //              Color col = sample_buffer[y * width + x];
        //              for (int k = 0; k < 3; ++k) {
        //                  this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        //              }
        //          }
        //      }
        // TODO: Task 2: You will likely want to update this function for supersampling support
        for (int x = 0; x < width; ++x) {
            for (int y = 0; y < height; ++y) {
                Color col(0, 0, 0);
                for (int k = 0; k < sample_rate; k++) {
                    col += sample_buffer[(y * width + x) * sample_rate + k];
                }
                col *= (1.0 / sample_rate);
                for (int k = 0; k < 3; ++k) {
                    this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
                }
            }
        }
    }

    Rasterizer::~Rasterizer() {}


}// CGL