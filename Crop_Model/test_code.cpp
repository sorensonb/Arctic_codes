/*
 *
 *
 *
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <codecvt>
#include <locale>
#include <string>
#include <math.h>

#include "test_header.h"

//double ConvertDegreesToRadians(double);
//int getMukeyfromLatLon(float, float);
//double find_distance(double, double, double, double);

int main(int argc, char *argv[])
{
    using std::cout;
    using std::wcout;
    using std::endl;

    double inlat;
    double inlon;

    //int f_mukey;
    std::wstring f_mukey;

    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    //
    // Prompt the user to input a lat and lon 
    //
    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

    if (argc < 3) 
    {
        cout << "SYNTAX: ./test_exec lat lon" << endl;
        return 0;
    }
    else
    {
        sscanf(argv[1], "%lf", &inlat);
        sscanf(argv[2], "%lf", &inlon);
        //inlat = atod(argv[1]);
        //inlon = atod(argv[2]);
    }

    f_mukey = getMukeyfromLatLon(inlat, inlon);

    wcout << L"Mukey = " << f_mukey << endl;    
    
    return 0;

}

double ConvertDegreesToRadians(double degree)
{
    double pi = 3.14159265359;
    return (degree * (pi / 180));
}

double find_distance(double rad_inlat, double rad_inlon, double rad_slat, double rad_slon)
{
    double calc_a;
    double calc_c;
    double dlat;
    double dlon;
    double distance;
    float earth_rad = 6371.;

    // Calculate the differences between the two input points
    // ------------------------------------------------------
    dlon = rad_slat - rad_inlat;
    dlat = rad_slon - rad_inlon;

    calc_a = pow(sin(dlat / 2), 2) + cos(rad_inlat) * cos(rad_slat) * pow(sin(dlon / 2), 2);
    calc_c = 2 * atan2(sqrt(calc_a) , sqrt(1 - calc_a));
    distance = earth_rad * calc_c;

    return distance;

}

std::wstring getMukeyfromLatLon(double inlat, double inlon) 
{
    using std::cout;
    using std::wcout;
    using std::endl;

    double rad_inlat;
    double rad_inlon;

    double slat;
    double slon;
    std::string s_mukey;
    std::wstring w_mukey;
    //int s_mukey;

    double rad_slat;
    double rad_slon;
    double distance;

    //int f_mukey;
    std::wstring f_mukey;
    double f_lat;
    double f_lon;
    double f_dist;
    int count;

    double max_distance;

    std::vector<std::string> row;
    std::string line, word, temp;

    // Convert the input latitude and longitude to radians 
    // (to save time later)
    // ----------------------------------------------------
    rad_inlat = ConvertDegreesToRadians(inlat);
    rad_inlon = ConvertDegreesToRadians(inlon);

    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    //
    // Determine the matching mukey for the input lat lon
    //
    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 
    // Open the file
    // -------------
    std::ifstream ifs("test_mukey.csv", std::ifstream::in);

    // Initialize the variables
    // ------------------------
    count = 0;
    f_lat = 89.99;
    f_lon = 179.99;
    f_dist = 99999.;
    s_mukey = "";
    //s_mukey = 0;
    while (!ifs.eof())
    {
        row.clear();
        
        getline(ifs, line);

        std::stringstream s(line);

        // Using the delimiter ',', get each word from the line
        // ----------------------------------------------------
        while (getline(s, word, ',')) 
        {
            row.push_back(word);
        }

        if(count != 0)
        {
            // Extract the latitude, longitude, and mukey from this line
            // ---------------------------------------------------------
            slat = stod(row[0]);
            slon = stod(row[1]);
            s_mukey = row[2];
            w_mukey = std::wstring_convert<std::codecvt_utf8<wchar_t>>().from_bytes(s_mukey);
            //s_mukey = stoi(row[2]);

            // Convert the grid lat/lon to radians
            // -----------------------------------
            rad_slat = ConvertDegreesToRadians(slat);
            rad_slon = ConvertDegreesToRadians(slon);

            // Determine the Great Circle distance between the points
            // ------------------------------------------------------
            distance = find_distance(rad_inlat, rad_inlon, rad_slat, rad_slon);

            // Determine if the distance from the current grid lat/lon to
            // the input lat/lon is less than the current "matching" point
            // -----------------------------------------------------------
            if(distance < f_dist) 
            {
                // Since this grid point is closer to the input lat/lon 
                // than the current matching point, set the current matching
                // point to this point.
                // --------------------------------------------------------- 
                f_lat = slat;
                f_lon = slon;
                f_dist = distance;
                f_mukey = w_mukey;
            }

        }

        count++;
    }

    // Check if the distance between the grid point and the input point is
    // too big
    // -------------------------------------------------------------------
    max_distance = 1.0;
    if(f_dist > max_distance) 
    {
        cout << "WARNING: Identified MUKEY point is far away from " << 
                "input point" << endl;
    }
    cout << "Input values: " << endl;
    cout << "   Lat =   " << inlat << endl;
    cout << "   Lon =   " << inlon << endl;
    cout << "Matching values: " << endl;
    cout << "   Lat =   " << f_lat << endl;
    cout << "   Lon =   " << f_lon << endl;
    cout << "   Dist =  " << f_dist << endl;
    wcout << "   Mukey = " << f_mukey << endl;

    //cout << "lines = " << count << endl;


    // Close the file
    // --------------
    ifs.close();

    return f_mukey;
}
