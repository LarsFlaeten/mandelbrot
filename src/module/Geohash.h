#ifndef _EXAMPLE_H_
#define _EXAMPLE_H_

#include <string>

namespace geohash {

    class Geohash {
		public:
			struct GeoBoxDimension {

   				double width;
    			double height;
 
			};
            struct GeoCoord {
    
                double latitude;
                double longitude;
    
                double north;
                double east;
                double south;
                double west;

    
            };

        public:
            Geohash();
            ~Geohash();

            
            // Encode a geohash
            // Accuracy radius is in meters
            // Higher radius (lower accuracy) means shorter hash
            std::string encode(double lat, double lon, double accRadius) const;

            // Decode a geohash
            // The size of the bounding box will be according to the accuracy of the hash
            GeoCoord    decode(const std::string& hash) const;


            int accuracyRadiusToInt(double accuracyRadius) const;
 
            GeoBoxDimension intToDegrees(int level) const;

    };
}

#endif
