#include <chrono>
#include <iostream>
#include <opencv2/opencv.hpp>

std::vector<cv::Point2f> control_points;

int beziernum = 4;
double step = 0.001;
void mouse_handler(int event, int x, int y, int flags, void *userdata) 
{
    if (event == cv::EVENT_LBUTTONDOWN && control_points.size() < beziernum)
    {
        std::cout << "Left button of the mouse is clicked - position (" << x << ", "
        << y << ")" << '\n';
        control_points.emplace_back(x, y);
    }     
}

void naive_bezier(const std::vector<cv::Point2f> &points, cv::Mat &window) 
{
    auto &p_0 = points[0];
    auto &p_1 = points[1];
    auto &p_2 = points[2];
    auto &p_3 = points[3];

    for (double t = 0.0; t <= 1.0; t += step)
    {
        auto point = std::pow(1 - t, 3) * p_0 + 3 * t * std::pow(1 - t, 2) * p_1 +
                 3 * std::pow(t, 2) * (1 - t) * p_2 + std::pow(t, 3) * p_3;

        window.at<cv::Vec3b>(point.y, point.x)[2] = 255;
    }
}

cv::Point2f recursive_bezier(const std::vector<cv::Point2f> &control_points, float t) 
{
    // TODO: Implement de Casteljau's algorithm
	int n = control_points.size();
	std::vector<cv::Point2f> res;
	float t1 = 1 - t;
	for (int i = 0; i < n - 1; i++)
	{
		res.push_back(t * control_points.at(i) + t1 * control_points.at(i + 1));
	}
	if (res.size() == 1)
	{
		return *res.begin();
	}
	else
	{
		return recursive_bezier(res, t);
	}

}

void bezier(const std::vector<cv::Point2f> &control_points, cv::Mat &window) 
{
    // TODO: Iterate through all t = 0 to t = 1 with small steps, and call de Casteljau's 
    // recursive Bezier algorithm.

	for (double t = 0.0; t <= 1.0; t += step)
	{
		cv::Point2f p = recursive_bezier(control_points, t);
		//window.at<cv::Vec3b>(p.y, p.x)[1] = 255;

		///*
		//点 p 最相邻的四个点
		cv::Point2i p0(p.x - std::floor(p.x) < 0.5 ? std::floor(p.x) : std::ceil(p.x),
			p.y - std::floor(p.y) < 0.5 ? std::floor(p.y) : std::ceil(p.y));
		std::vector<cv::Point2i> ps = { p0, cv::Point2i(p0.x - 1, p0.y),
			cv::Point2i(p0.x, p0.y - 1), cv::Point2i(p0.x - 1, p0.y - 1),
		};

		// 点 p 到相邻四个点的中心点的距离
		float sum_d = 0.f;
		float max_d = sqrt(2);
		std::vector<float> ds = {};
		for (int i = 0; i < 4; i++) {
			cv::Point2f cp(ps[i].x + 0.5f, ps[i].y + 0.5f);
			float d = max_d - std::sqrt(std::pow(p.x - cp.x, 2) + std::pow(p.y - cp.y, 2));
			ds.push_back(d);
			sum_d += d;
		};

		// assign colors
		for (int i = 0; i < 4; i++) {
			float k = ds[i] / sum_d;
			window.at<cv::Vec3b>(ps[i].y, ps[i].x)[1] = std::min(255.f, window.at<cv::Vec3b>(ps[i].y, ps[i].x)[1] + 255.f * k);
		};
		//*/
		
	}

}

int main() 
{
    cv::Mat window = cv::Mat(700, 700, CV_8UC3, cv::Scalar(0));
    cv::cvtColor(window, window, cv::COLOR_BGR2RGB);
    cv::namedWindow("Bezier Curve", cv::WINDOW_AUTOSIZE);

    cv::setMouseCallback("Bezier Curve", mouse_handler, nullptr);

    int key = -1;
    while (key != 27) 
    {
        for (auto &point : control_points) 
        {
            cv::circle(window, point, 3, {255, 255, 255}, 3);
        }

        if (control_points.size() == beziernum)
        {
            //naive_bezier(control_points, window);
            bezier(control_points, window);

            cv::imshow("Bezier Curve", window);
            cv::imwrite("my_bezier_curve.png", window);
            //key = cv::waitKey(20);

            //return 0;
        }

        cv::imshow("Bezier Curve", window);
        key = cv::waitKey(20);
		if (key == -1)
		{

		}
		else if (key == 0)//F1
		{
			//清空

			//window.release();
			window.setTo(cv::Scalar(0));
			control_points.clear();
			/*w*indow = cv::Mat(700, 700, CV_8UC3, cv::Scalar(0));
			cv::cvtColor(window, window, cv::COLOR_BGR2RGB);
			cv::namedWindow("Bezier Curve", cv::WINDOW_AUTOSIZE);

			cv::setMouseCallback("Bezier Curve", mouse_handler, nullptr);*/
		}
		else
		{
			if (key == 13)//Enter键
			{
				control_points.clear();
				window.setTo(cv::Scalar(0));

				std::cout << "请输入贝塞尔曲线的控制点"<<std::endl;
				std::cin >> beziernum ;
			}
			std::cout << key << std::endl;
		}
    }

return 0;
}
