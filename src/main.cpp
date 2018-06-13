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


// initializing lane 
int lane = 1;
//initializing reference speed in mph
double ref_speed = 0; // max should be 49.5

  h.onMessage([&max_s,&lane, &ref_speed, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
		
		double pos_x;
		double pos_y;
		double angle;
		int path_size = previous_path_x.size();
		// debug line cout << "debugging length of path size " << path_size << endl;
		// initializing vector of pts and reference values

		vector<double> ptsx;
		vector<double> ptsy;
		double ref_x = car_x;
		double ref_y = car_y;
		double ref_yaw = deg2rad(car_yaw);		
		
		

		if (path_size>0){
			car_s = end_path_s;
			car_d = end_path_d;
		}

		// lane calc
		if (car_d > 0 && car_d < 4) { lane = 0;}
		else if (car_d > 4 && car_d < 8) { lane = 1;}	
		else if (car_d > 8 && car_d < 12) { lane = 2;}
		
		cout << "car lane is " << car_d << " lane = " << lane << endl;		
		int current_lane = 2+4*lane;


		// indicator for vehicle too close to ego car
		bool too_close = false;
		bool slow_down = false;

		// method for checking lanes and iterating over 
		vector<bool> lane_check = { false, false, false};
		// indicator for checking safe transitionary state change
		vector<bool> transition_check = { false, false, false};
		for (int i = 0; i<sensor_fusion.size(); i++) // looping for lane safety
		{
			double vx = sensor_fusion[i][3];
			double vy = sensor_fusion[i][4];
			double check_speed = sqrt(vx*vx+vy*vy);
			double check_car_s = sensor_fusion[i][5];
			double other_d = sensor_fusion[i][6];

			//previous path size used to extrapolate s from speed
			check_car_s += ((double)path_size*.02*check_speed);
			
			// checking if cars prompt lane change or slowdown
			if( (check_car_s > car_s) && ((check_car_s-car_s) < 50))
			{	
				for(int l = 0; l<lane_check.size(); l++) 
				{			
					double lane_compare = 2.0+4.0*(double)l;
					if((other_d>lane_compare-2) && (other_d<lane_compare+2)) {lane_check[l]=true;} 	
										
				}
			} // end if lane_check 

			// checking if nearby lanes are possible transitonary states
			if( (check_car_s > car_s-8) && ((check_car_s-car_s) < 10))
			{	
				for(int l = 0; l<transition_check.size(); l++) 
				{			
					double lane_compare = 2.0+4.0*(double)l;
					if((other_d>lane_compare-2) && (other_d<lane_compare+2)) {transition_check[l]=true;} 	
										
				}
			} // end if transition_check

			} // end for loop


		cout << " lane check :" << lane_check[0] << lane_check[1] << lane_check[2] << endl;
		cout << " transition check :" << transition_check[0] << transition_check[1] << transition_check[2] << endl;
		if (lane_check[lane]) {too_close = true;}
		//else {too_close=false;}
		
		int possible_lane = 3;
		if (too_close) {
			possible_lane =3;	
			for (int l = 0; l<lane_check.size(); l++) 
			{
				if (!lane_check[l]) {possible_lane = l;}
			}
			
			if ((lane % 2 == 0) && (!lane_check[1])) 
			{
				lane = 1;
				current_lane = 2+4*lane; 
				cout << "changing to lane " << lane << endl;
			}
			if ((lane == 0) && (lane_check[1]))
			{   
	 			slow_down = true;
				if ((!transition_check[1])&&(!lane_check[2]))
				{
				lane = 1;
				current_lane = 2+4*lane; 
				cout << "changing to lane " << lane << endl;
				}	
			}
			else if ((lane == 2) && (lane_check[1]))
			{   
	 			slow_down = true;
				if ((!transition_check[1])&&(!lane_check[0]))
				{
				lane = 1;
				current_lane = 2+4*lane; 
				cout << "changing to lane " << lane << endl;
				}				
			}
			else if (lane==1 && possible_lane !=3) // add criteria to slow down in failed if statements
			{ 				
				lane = possible_lane;
				current_lane = 2+4*lane; 
				cout << "changing to lane " << lane << endl;
			}

			else if (possible_lane==3) {slow_down = true;} 		
		
		}
		
		if(slow_down) 
		{
			ref_speed -= .224;
			cout << " slowing down" <<endl;
		}
				
		else if (ref_speed<49.0) 
		{
		    ref_speed += .224; 
			cout << "Speeding up, ref speed: "<< ref_speed << endl;
		}

		if(ref_speed>49.0) {
		ref_speed = 49.0; 
		cout << "Maintaining speed " << endl;
		}

		
			//cout << "ref speed is " << ref_speed << endl;
		//cout << " too close flag is " << too_close << endl;

		//int randint = rand() % 2;		
		//if (lane_change) { lane = randint*2;}

		// if pathsize is too small, use reference values to make path tangent to car
		if (path_size < 2) {
		double prev_car_x = car_x - cos(car_yaw);		
		double prev_car_y = car_y - sin(car_yaw);
		
		ptsx.push_back(prev_car_x);
		ptsx.push_back(car_x);

		ptsy.push_back(prev_car_y);
		ptsy.push_back(car_y);
				
		}
		
		else {
		
		// else use previous path end point
		ref_x = previous_path_x[path_size-1];
		ref_y = previous_path_y[path_size-1];
		double ref_prev_x = previous_path_x[path_size-2];
		double ref_prev_y = previous_path_y[path_size-2];
		
		ref_yaw = atan2(ref_y-ref_prev_y, ref_x-ref_prev_x);				
		
		ptsx.push_back(ref_prev_x);
		ptsy.push_back(ref_prev_y);
		
		ptsx.push_back(ref_x);
		ptsy.push_back(ref_y);

		}
		 // debug line cout << "debugging length of ptsx " << ptsx.size() << endl;

		// using frenet coordinates to space subsequent pts 30 m ahead

		
		vector<double> next_wps0 = getXY(car_s+30,current_lane, map_waypoints_s,map_waypoints_x, map_waypoints_y);		
		vector<double> next_wps1 = getXY(car_s+60,current_lane, map_waypoints_s,map_waypoints_x, map_waypoints_y);
		vector<double> next_wps2 = getXY(car_s+90,current_lane, map_waypoints_s,map_waypoints_x, map_waypoints_y);

		ptsx.push_back(next_wps0[0]);
		ptsx.push_back(next_wps1[0]);
		ptsx.push_back(next_wps2[0]);

		ptsy.push_back(next_wps0[1]);
		ptsy.push_back(next_wps1[1]);
		ptsy.push_back(next_wps2[1]);
		
		// debug line cout << "debugging length of ptsx " << ptsx.size() << endl; //debug size

		// shifting car's reference angle to 0
		for (int i = 0; i<ptsx.size(); i++){

		double shift_x = ptsx[i]-ref_x ;
		double shift_y = ptsy[i]-ref_y ;
	
		ptsx[i] = (shift_x*cos(0-ref_yaw) - shift_y*sin(0-ref_yaw)) ; 
		ptsy[i] = (shift_x*sin(0-ref_yaw) + shift_y*cos(0-ref_yaw)) ;		
		}

		// initializing spline and setting polynomial from ptsx,ptsy
		tk::spline s;
		s.set_points(ptsx,ptsy);
		
		vector<double> next_x_vals;
          	vector<double> next_y_vals;

		
		// start by setting previous path pts to next x,y planned coordinates
		for (int i = 0; i<previous_path_x.size();i++){
		
		next_x_vals.push_back(previous_path_x[i]);
		next_y_vals.push_back(previous_path_y[i]);
		
		} // using previous path references 

		// debug line cout << "debugging length of next_x_vals " << next_x_vals.size() << endl; //debug line

		// calculate spacing of spline points
		double target_x = 30.0;
		double target_y = s(target_x);
		double target_distance = sqrt((target_x*target_x)+(target_y*target_y));

		double x_addon = 0.0;

		// using equation N * latency_interval * speed = distance
		// N = distance/(speed*latency_interval
		double N = (target_distance / (.02*ref_speed/2.24));
		// debug line cout << "N value is " << N << endl;


 		// filling out remaining points from spline
		for (int i = 0; i < 50-previous_path_x.size(); i++) {


		double x_point = x_addon + (target_x/N);
		double y_point = s(x_point);
		// debug line cout << "x value is " << x_point << endl;
		// debug line cout << "y value is " << y_point << endl;		
		x_addon = x_point;

		
		double x_ref = x_point;
		double y_ref = y_point;

		// rotate back from reference angle 0, reverting shift above.
		x_point = (x_ref*cos(ref_yaw) - y_ref*sin(ref_yaw)); 
		y_point = (x_ref*sin(ref_yaw) + y_ref*cos(ref_yaw)); 	


		x_point += ref_x;   // different ref_x and x_ref
		y_point += ref_y;   // different ref_y and y_ref


		next_x_vals.push_back(x_point);
		next_y_vals.push_back(y_point);
		
		}
		// debug line cout << "debugging length of next_x_vals " << next_x_vals.size() << endl; //debug line


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
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
