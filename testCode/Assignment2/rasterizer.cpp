// clang-format off
//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>


rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f> &cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}

static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f* v)
{
	float c1 = (x*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*y + v[1].x()*v[2].y() - v[2].x()*v[1].y()) / (v[0].x()*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*v[0].y() + v[1].x()*v[2].y() - v[2].x()*v[1].y());
	float c2 = (x*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*y + v[2].x()*v[0].y() - v[0].x()*v[2].y()) / (v[1].x()*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*v[1].y() + v[2].x()*v[0].y() - v[0].x()*v[2].y());
	float c3 = (x*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*y + v[0].x()*v[1].y() - v[1].x()*v[0].y()) / (v[2].x()*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*v[2].y() + v[0].x()*v[1].y() - v[1].x()*v[0].y());
	return { c1,c2,c3 };
}

static bool insideTriangle(float x, float y, const Vector3f* _v)
{   
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
	//测试点是否在三角形内
	auto[alpha, beta, gamma] = computeBarycentric2D(x, y, _v);
	return (alpha >= 0 && alpha <= 1 &&
		beta >= 0 && beta <= 1 &&
		gamma >= 0 && gamma <= 1);

}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type, bool ssaa)
{
    auto& buf = pos_buf[pos_buffer.pos_id];
    auto& ind = ind_buf[ind_buffer.ind_id];
    auto& col = col_buf[col_buffer.col_id];

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (auto& i : ind)
    {
        Triangle t;
        Eigen::Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f)
        };
        //Homogeneous division
        for (auto& vec : v) {
            vec /= vec.w();
        }
        //Viewport transformation
        for (auto & vert : v)
        {
            vert.x() = 0.5*width*(vert.x()+1.0);
            vert.y() = 0.5*height*(vert.y()+1.0);
            vert.z() = vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
        }

        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];

        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);
		if (!ssaa)
		{
			rasterize_triangle(t);
		}
		else
		{
			rasterize_triangle_ssaa(t);
		}
    }
	if (ssaa)
	{
		updateframeBuffer_ssaa();
	}
}

//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle& t) 
{
    auto v = t.toVector4();
    
    // TODO : Find out the bounding box of current triangle.
    // iterate through the pixel and find if the current pixel is inside the triangle

    // If so, use the following code to get the interpolated z value.
    //auto[alpha, beta, gamma] = computeBarycentric2D(x, y, t.v);
    //float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
    //float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
    //z_interpolated *= w_reciprocal;

    // TODO : set the current pixel (use the set_pixel function) to the color of the triangle (use getColor function) if it should be painted.
	//找出三角形的边界
	float xmin, ymin, xmax, ymax;
	xmin = xmax = t.v[0].x();
	ymin = ymax = t.v[0].y();

	xmin = std::min(xmin, t.v[1].x()); xmax = std::max(xmax, t.v[1].x());
	ymin = std::min(ymin, t.v[1].y()); ymax = std::max(ymax, t.v[1].y());
	xmin = std::min(xmin, t.v[2].x()); xmax = std::max(xmax, t.v[2].x());
	ymin = std::min(ymin, t.v[2].y()); ymax = std::max(ymax, t.v[2].y());

	for (int i = xmin; i < xmax; i++)
	{
		for (int j = ymin; j < ymax; j++)
		{
			if (insideTriangle(i + 0.5, j + 0.5, t.v))//判断是否在三角形内
			{
				auto[alpha, beta, gamma] = computeBarycentric2D(i + 0.5, j + 0.5, t.v);
				float w_reciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
				float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
				z_interpolated *= w_reciprocal;
				auto ind = get_index(i, j);
				if (z_interpolated < depth_buf[ind])
				{
					depth_buf[ind] = z_interpolated;
					Eigen::Vector3f point(i, j, 0);
					set_pixel(point, t.getColor());
				}
			}

		}
	}
}

void rst::rasterizer::rasterize_triangle_ssaa(const Triangle & t)
{
	auto v = t.toVector4();

	//找出三角形的边界
	float xmin, ymin, xmax, ymax;
	xmin = xmax = t.v[0].x();
	ymin = ymax = t.v[0].y();

	xmin = std::min(xmin, t.v[1].x()); xmax = std::max(xmax, t.v[1].x());
	ymin = std::min(ymin, t.v[1].y()); ymax = std::max(ymax, t.v[1].y());
	xmin = std::min(xmin, t.v[2].x()); xmax = std::max(xmax, t.v[2].x());
	ymin = std::min(ymin, t.v[2].y()); ymax = std::max(ymax, t.v[2].y());

	for (int i = xmin; i < xmax; i++)
	{
		for (int j = ymin; j < ymax; j++)
		{
			auto ind = get_index(i, j);
			for (int x = 0; x < ssaa_w; x += 1)
			{
				for (int y = 0; y < ssaa_h; y += 1)
				{
					if (insideTriangle(i + x * ssaa_pixel_step_w + ssaa_pixel_step_w / 2.0, j + y * ssaa_pixel_step_h + ssaa_pixel_step_h / 2.0, t.v))//判断是否在三角形内
					{
						auto[alpha, beta, gamma] = computeBarycentric2D(i + x * ssaa_pixel_step_w + ssaa_pixel_step_w / 2.0, j + y * ssaa_pixel_step_h + ssaa_pixel_step_h / 2.0, t.v);
						float w_reciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
						float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
						z_interpolated *= w_reciprocal;
						if (z_interpolated < depth_buf_ssaa[ind][x * ssaa_w + y])
						{
							depth_buf_ssaa[ind][x * ssaa_w + y] = z_interpolated;
							frame_buf_ssaa[ind][x * ssaa_w + y] = t.getColor();
						}
					}
				}
			}		

		}
	}
	//for (int x = xmin; x < xmax; x++)
	//{
	//	for (int y = ymin; y < ymax; y++) {

	//		Eigen::Vector3f point(x, y, 1.0);
	//		//child pixel
	//		int inside_count = 0;
	//		int update_depth = 0;
	//		int index = 0;
	//		for (float i = 0.25; i < 1.0; i += 0.5) {
	//			for (float j = 0.25; j < 1.0; j += 0.5)
	//			{
	//				if (insideTriangle(x + i, y + j, t.v)) {
	//					auto[alpha, beta, gamma] = computeBarycentric2D(x + i, y + j, t.v);
	//					float w_reciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
	//					float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
	//					z_interpolated *= w_reciprocal;
	//					if (z_interpolated < depth_buf_ssaa[get_index(x, y)][index]) {

	//						point << x + i, y + j, 1.0;
	//						//set_pixel_ssaa(point, index, t.getColor());
	//						frame_buf_ssaa[get_index(x, y)][index] = t.getColor();
	//						depth_buf_ssaa[get_index(x, y)][index] = z_interpolated;
	//						inside_count++;
	//						update_depth += z_interpolated;
	//					}

	//				}
	//				index++;
	//			}

	//		}

	//	}
	//}
}

void rst::rasterizer::updateframeBuffer_ssaa()
{
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			Eigen::Vector3f color = { 0., 0., 0. };
			int index = get_index(i, j);
			bool isEdge{ false };
			Eigen::Vector3f temp{ frame_buf_ssaa[index][0] };
			for (int x = 0; x < ssaa_w; x += 1)
			{
				for (int y = 0; y < ssaa_h; y += 1)
				{
					if (temp != frame_buf_ssaa[index][x * ssaa_w + y])
					{
						isEdge = true;
					}
					color += frame_buf_ssaa[index][x * ssaa_w + y];
				}
			}
			Vector3f point(i, j, 0);
			set_pixel(point, color / (ssaa_w * ssaa_h));
			if (isEdge)
			{
				//set_pixel(point, Vector3f(255,255,255));
			}
	
		}
	}
	//for (int x = 0; x < width; x++) {
	//	for (int y = 0; y < height; y++)
	//	{
	//		Eigen::Vector3f color(0, 0, 0);
	//		for (int i = 0; i < 4; i++)
	//			color += frame_buf_ssaa[get_index(x, y)][i];
	//		color /= 4;
	//		set_pixel(Eigen::Vector3f(x, y, 1.0f), color);
	//	}

	//}
}

void rst::rasterizer::set_model(const Eigen::Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
		for(int i = 0; i < frame_buf_ssaa.size(); i++) {
			std::fill(frame_buf_ssaa[i].begin(), frame_buf_ssaa[i].end(), Eigen::Vector3f{ 0, 0, 0 });
		}
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
		for (int i = 0; i < depth_buf_ssaa.size(); i++) {
			std::fill(depth_buf_ssaa[i].begin(), depth_buf_ssaa[i].end(), std::numeric_limits<float>::infinity());
		}
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);
	frame_buf_ssaa.resize(w * h); 
	depth_buf_ssaa.resize(w * h);
	for (int i = 0; i < frame_buf_ssaa.size(); i++) {
		frame_buf_ssaa[i].resize(ssaa_w * ssaa_h);
	}
	for (int i = 0; i < depth_buf_ssaa.size(); i++) {
		depth_buf_ssaa[i].resize(ssaa_w * ssaa_h);
	}
}

int rst::rasterizer::get_index(int x, int y)
{
    return (height-1-y)*width + x;
}

void rst::rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    //old index: auto ind = point.y() + point.x() * width;
    auto ind = (height-1-point.y())*width + point.x();
    frame_buf[ind] = color;

}

// clang-format on