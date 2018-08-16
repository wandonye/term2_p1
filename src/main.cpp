#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

bool checkCloseFront(double check_car_speed, double check_car_s, 
						double prev_size, double car_speed, double car_s, double front_margin) {
	check_car_s+=prev_size*0.02*check_car_speed;
	if ((check_car_s>car_s)&&((check_car_s-car_s)<front_margin*car_speed/50)) {
		return true;
	}
	return false;
}

bool checkSideCarSafe(double check_car_speed, double check_car_s, 
						double prev_size, double car_speed, double car_s, double front_margin, double back_margin) {
	check_car_s+=prev_size*0.02*check_car_speed;
  double delta = check_car_s - car_s;
  if ((delta<front_margin*car_speed/50) && (delta>-back_margin*check_car_speed/(car_speed+1))) {
    // the car being checked is within a window of the current car
    // the boundary of the window in the back depends on the relative speed of the two car
    return false;
  }
	return true;
}

// All the possible states of the car
enum carStates
{
    keepLane,
    willChangeLeft,
    willChangeRight,
    isChangingLeft,
    isChangingRight
};

// initiate car state as lane keeping
carStates cur_state = carStates::keepLane;

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  int lane = 1;
  double ref_speed = 0.0;

  h.onMessage([&lane,&ref_speed,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          int prev_size = previous_path_x.size();
          // reuse previous plan
          if(prev_size > 0) {
            car_s = end_path_s;
          }

          // if there is a close car in front in the same lane
          bool too_close_lane = false;
          // if there is a close car in front in the left lane
          bool is_ok_to_left = (car_d>4);
          // if there is a close car in front in the right lane
          bool is_ok_to_right = (car_d<8);
                
          // distance threshold for car in the same lane
          double front_margin = 25.0;
          // distance threshold for car in other lanes
          double lane_change_front_margin = 25.0;
          double lane_change_back_margin = 10.0;

          // speed of the car within threshold (if any)
          double front_car_speed = 100.0;
          double left_car_speed = 100.0;
          double right_car_speed = 100.0;

          // sensor fusion
          for (int i=0; i<sensor_fusion.size(); i++) {   // Evaluate each car
            float d = sensor_fusion[i][6];
            double vx = sensor_fusion[i][3];
            double vy = sensor_fusion[i][4];
            // car is in lane number=lane
            double speed = sqrt(vx*vx+vy*vy);

            if (d<(4*lane)) {
              // sensed car is in the left lane
              if (lane!=0) {
                // no need to check left lane if it's already the leftmost
                // cout<<"left lane:"<<sensor_fusion[i]<<endl;
                if (!checkSideCarSafe(speed,sensor_fusion[i][5],
                                prev_size, car_speed, car_s, lane_change_front_margin, lane_change_back_margin)) {
                  left_car_speed = min(speed, left_car_speed);
                  is_ok_to_left = false;
                }
              }
            } else if (d<(4*lane+4)) {
              // cout<<"current lane:"<<sensor_fusion[i]<<endl;
              if (checkCloseFront(speed, sensor_fusion[i][5],
                              prev_size, car_speed, car_s, front_margin)) {
                front_car_speed = min(speed, front_car_speed);
                too_close_lane = true;
              }
            } else if (d<(4*lane+8)) {
              // cout<<"right lane:"<<sensor_fusion[i]<<endl;
              if (!checkSideCarSafe(speed, sensor_fusion[i][5],
                              prev_size, car_speed, car_s, lane_change_front_margin, lane_change_back_margin)) {
                right_car_speed = min(speed, right_car_speed);
                is_ok_to_right = false;
              }
            }
          }
          // cout<<"end of one scan"<<endl;
          
          // path anchor points
          vector<double> ptsx;
          vector<double> ptsy;

          // ref_x, ref_y: the starting point of the new calculation
          double ref_x = car_x;
          double ref_y = car_y;
          double rot = deg2rad(car_yaw);
          if (prev_size<2) {
            // when less the 2 points left
            // use car's current location and one step backward as 
            // two initial points for this round of calculation
            double prev_car_x = car_x-cos(rot);
            double prev_car_y = car_y-sin(rot);

            ptsx.push_back(prev_car_x);
            ptsx.push_back(car_x);
            ptsy.push_back(prev_car_y);
            ptsy.push_back(car_y);
          } else {
            // otherwise, use last two points of the previous plan 
            // as initial points
            ref_x = previous_path_x[prev_size-1];
            ref_y = previous_path_y[prev_size-1];

            double ref_x_prev = previous_path_x[prev_size-2];
            double ref_y_prev = previous_path_y[prev_size-2];
            rot = atan2(ref_y-ref_y_prev,ref_x-ref_x_prev);

            ptsx.push_back(ref_x_prev);
            ptsx.push_back(ref_x);
            ptsy.push_back(ref_y_prev);
            ptsy.push_back(ref_y);
          }

          // State transition
          switch (cur_state)
          {
            case carStates::keepLane:
              // cout<<"keepLane"<<endl;
              if (too_close_lane) {
                if (is_ok_to_left&&is_ok_to_right) {
                  // if both lane are available, choose the faster lane
                  cur_state = (right_car_speed>left_car_speed)?carStates::willChangeRight:carStates::willChangeLeft;
                } else if (is_ok_to_left) {
                  // prepare left lane change
                  cur_state = carStates::willChangeLeft;
                } else if (is_ok_to_right) {
                  // prepare right lane change
                  cur_state = carStates::willChangeRight;
                } else {
                  if (ref_speed>front_car_speed) { ref_speed -= 0.224; }
                }
              } else if (ref_speed<49.75) {
                ref_speed += 0.224;
              }
              break;
            case carStates::willChangeLeft:
              // cout<<"willChangeLeft"<<endl;
              lane = lane - 1;
              cur_state = carStates::isChangingLeft;
              break;
            case carStates::willChangeRight:
              // cout<<"willChangeRight"<<endl;
              lane = lane + 1;
              cur_state = carStates::isChangingRight;
              break;
            case carStates::isChangingLeft:
              // cout<<"isChangingLeft"<<endl;
              if (car_d<(4*lane+4)&&car_d>4*lane) {
                cur_state = carStates::keepLane;
              }
              break;
            case carStates::isChangingRight:
              // cout<<"isChangingRight"<<endl;
              if (car_d<(4*lane+4)&&car_d>4*lane) {
                cur_state = carStates::keepLane;
              }
              break;
            default:
              break;
          }

          vector<double> next_wp0 = getXY(car_s+30,2+4.0*lane,map_waypoints_s,map_waypoints_x,map_waypoints_y);
          vector<double> next_wp1 = getXY(car_s+60,2+4.0*lane,map_waypoints_s,map_waypoints_x,map_waypoints_y);
          vector<double> next_wp2 = getXY(car_s+90,2+4.0*lane,map_waypoints_s,map_waypoints_x,map_waypoints_y);

          ptsx.push_back(next_wp0[0]);
          ptsx.push_back(next_wp1[0]);
          ptsx.push_back(next_wp2[0]);

          ptsy.push_back(next_wp0[1]);
          ptsy.push_back(next_wp1[1]);
          ptsy.push_back(next_wp2[1]);

          for (int i=0; i<ptsx.size(); i++) {
            double dx = ptsx[i]-ref_x;
            double dy = ptsy[i]-ref_y;
            
            ptsx[i] = dx*cos(0-rot)-dy*sin(0-rot);
            ptsy[i] = dx*sin(0-rot)+dy*cos(0-rot);
          }
          
          // create a spline
          tk::spline s;
          // set anchor points of the spline
          s.set_points(ptsx,ptsy);
          // now think of s as the function for the path curve, 
          // input x-coord, output y-coord
          // shifted ref point to the origin, i.e. input 0, output 0

          // use path points from previous iteration
          for (int i = 0; i < previous_path_x.size(); i++)
          {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          double target_x = 30.0;
          double target_y = s(target_x);

          // Euclidean distance
          double target_dist = sqrt(target_x*target_x+target_y*target_y);
          double x_add_on = 0;

          for(int i = 1; i < 50-prev_size; i++)
          {
            double N = target_dist/(0.02*ref_speed/2.24);
            double x_pt = x_add_on+target_x/N;
            double y_pt = s(x_pt);
            x_add_on = x_pt;

            double x_bak = x_pt;
            double y_bak = y_pt;

            x_pt = x_bak*cos(rot)-y_bak*sin(rot);
            y_pt = x_bak*sin(rot)+y_bak*cos(rot);

            x_pt += ref_x;
            y_pt += ref_y;

            next_x_vals.push_back(x_pt);
            next_y_vals.push_back(y_pt);
          }

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
