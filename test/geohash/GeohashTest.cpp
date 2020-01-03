#include <gtest/gtest.h>

#include <geohash/Geohash.h>
#include <cmath>

class GeohashTest : public ::testing::Test {

protected:
    GeohashTest() { ; }

    virtual ~GeohashTest() { ; }

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp() { ; }

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown() { ; }

	
	geohash::Geohash hasher;

};



TEST_F(GeohashTest, TestLevelToDimension)
{

	for(int i = 1; i < 13; ++i) {
		auto box = hasher.intToDegrees(i);
		

		auto level = hasher.accuracyRadiusToInt(box.height*110570);
	    ASSERT_EQ(i, level);	
	    //std::cout << i << ": " << box.width << " x " << box.height << ", " << box.height*110570 << "m" << ", Calculated level; " << hasher.accuracyRadiusToInt(box.height*110570) << std::endl;
		
	}

    // check tath more accurate than level 12 is cropped to 12:
    ASSERT_EQ(hasher.accuracyRadiusToInt(hasher.intToDegrees(13).height*110570), 12);

    // Check that < 1 is cropped to 1:
    ASSERT_EQ(hasher.accuracyRadiusToInt(hasher.intToDegrees(0).height*110570), 1);

} 

TEST_F(GeohashTest, TestRandomGeoHashOrg)
{

    auto hash = hasher.encode(10.89, 67.45, 100);
    ASSERT_EQ(hash, "t3zvyy7");
    
    hash = hasher.encode(0.001, 0.001, 10);
    ASSERT_EQ(hash, "s000000t");

    hash = hasher.encode(-25.382708, -49.265506, 0.01);
    ASSERT_EQ(hash, "6gkzwgjzn820");

    ASSERT_EQ(hasher.encode(28.5934, 56.5316, 1), "tm04wmm09");
    ASSERT_EQ(hasher.encode(47.5269, 94.2939, 1), "y07n3sn8r" );
    ASSERT_EQ(hasher.encode(-67.2061, -134.538, 1),  "1h03dt82j");
    ASSERT_EQ(hasher.encode(-19.5894, -39.5679, 1),  "7hebe9ge3");
    ASSERT_EQ(hasher.encode(4.2189, 7.91726, 0.1), "s0v8h0j0veb");

    auto geo = hasher.decode("s0v8h0j0veb");
    ASSERT_LT(fabs(geo.latitude-4.2189), 0.0001);//, 7.91726, 0.1), "s0v8h0j0veb");

   

} 

