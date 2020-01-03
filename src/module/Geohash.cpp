#include "Geohash.h"
#include <cmath>
namespace geohash {
    /* Normal 32 characer map used for geohashing */
    static std::string char_map =  "0123456789bcdefghjkmnpqrstuvwxyz";


    struct Interval {
    
        double high;
        double low;
    
    };

    Geohash::Geohash() {;}

    Geohash::~Geohash() {;}

    std::string Geohash::encode(double lat, double lng, double accuracyRadius) const {
        //std::string hash = "";
        
        int precision = accuracyRadiusToInt(accuracyRadius);
   
        char hash[precision];

        if(lat <= 90.0 && lat >= -90.0 && lng <= 180.0 && lng >= -180.0) {
        
            //hash = (char*)malloc(sizeof(char) * (precision + 1));
            hash[precision] = '\0';
        
            precision *= 5.0;
        
            Interval lat_interval = {90.0, -90.0};
            Interval lng_interval = {180.0, -180.0};

            Interval *interval;
            double coord, mid;
            int is_even = 1;
            unsigned int hashChar = 0;
            int i;
            for(i = 1; i <= precision; i++) {
         
                if(is_even) {
                
                    interval = &lng_interval;
                    coord = lng;                
                    
                } else {
                    
                    interval = &lat_interval;
                    coord = lat;   
                }
                
                mid = (interval->low + interval->high) / 2.0;
                hashChar = hashChar << 1;
                
                if(coord > mid) {
                    
                    interval->low = mid;
                    hashChar |= 0x01;
                    
                } else
                    interval->high = mid;
                
                if(!(i % 5)) {
                    
                    hash[(i - 1) / 5] = char_map[hashChar];
                    hashChar = 0;

                }
                
                is_even = !is_even;
            }
         
            
        }
                 
        

        return hash;
    }

    unsigned int index_for_char(char c, const std::string& string) {
    
        int index = -1;
        int string_amount = string.length();
        int i;
        for(i = 0; i < string_amount; i++) {
        
            if(c == string[i]) {
            
                index = i; 
                break;
            }
        
        }
    
        return index;
    }


    Geohash::GeoCoord    Geohash::decode(const std::string& hash) const {
        
        GeoCoord coordinate = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        int char_amount = hash.length();
        
        if(char_amount) {
            
            unsigned int char_mapIndex;
            Interval lat_interval = {90.0, -90.0};
            Interval lng_interval = {180.0, -180.0};
            Interval *interval;
        
            int is_even = 1;
            double delta;
            int i, j;
            for(i = 0; i < char_amount; i++) {
            
                char_mapIndex = index_for_char(hash[i], char_map);
                
                if(char_mapIndex < 0)
                    break;
            
                // Interpret the last 5 bits of the integer
                for(j = 0; j < 5; j++) {
                
                    interval = is_even ? &lng_interval : &lat_interval;
                
                    delta = (interval->high - interval->low) / 2.0;
                
                    if((char_mapIndex << j) & 0x0010)
                        interval->low += delta;
                    else
                        interval->high -= delta;
                
                    is_even = !is_even;
                }
            
            }
            
            coordinate.latitude = lat_interval.high - ((lat_interval.high - lat_interval.low) / 2.0);
            coordinate.longitude = lng_interval.high - ((lng_interval.high - lng_interval.low) / 2.0);
            
            coordinate.north = lat_interval.high;
            coordinate.east = lng_interval.high;
            coordinate.south = lat_interval.low;
            coordinate.west = lng_interval.low;
        }
         

        return coordinate;
    }


	int Geohash::accuracyRadiusToInt(double accuracyRadius) const {
		// https://stackoverflow.com/questions/13448595/geohash-string-length-and-accuracy
		// Bottom answer
		int level = std::floor(std::log2(5000000/accuracyRadius) / 2.5 + 1 );
		if(level < 1)
			return 1;
		else if (level > 12)
			return 12;
		else
			return level;


	}
/*
    int Geohash::accuracyRadiusToInt(double accuracyRadius) const {
		double deg = accuracyRadius/110570.0; // Sice of the radius in latitude
		
		int cuts = 0;
		while(deg < 180.0) {
			deg = deg * 2.0;
			cuts++;
		}

		int level = 2 * cuts / 5 + (cuts % 2 ? 0 : 1);
		
		return level;


    }
*/            
    Geohash::GeoBoxDimension Geohash::intToDegrees(int level) const {
       GeoBoxDimension dimensions = {360.0, 180.0};

    	if(level > 0) {

        	int lat_times_to_cut = level * 5 / 2;
        	int lng_times_to_cut = level * 5 / 2 + (level % 2 ? 1 : 0);

        	double width = 360.0;
        	double height = 180.0;

        	int i;
        	for(i = 0; i < lat_times_to_cut; i++)
            	height /= 2.0;

        	for(i = 0; i < lng_times_to_cut; i++)
            	width /= 2.0;

        	dimensions.width = width;
        	dimensions.height = height;

    	}

    	return dimensions; 
    }


}
