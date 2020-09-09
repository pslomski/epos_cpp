#include <catch.hpp>
#include <epos.hpp>

using namespace epos;
using namespace std;

TEST_CASE("Channel::getT")
{
	REQUIRE(Channel::getT(0) == 0.5);
	REQUIRE(Channel::getT(epos::NbrSAMP - 1) == 800);
}

TEST_CASE("Channel::getIdx")
{
	REQUIRE(Channel::getIdx(0.5) == 0);
	REQUIRE(Channel::getIdx(800) == (epos::NbrSAMP - 1));
}

TEST_CASE("Channel::init")
{
	Channel channel;
	REQUIRE(channel.data.size() == 0);

	channel.add(1.0);
	REQUIRE(channel.data.size() == 1);
	REQUIRE(channel.data[0] == 1.0);

	// clears data and only reserves storage of epos::NbrSAMP elems
	channel.init(epos::NbrSAMP);
	REQUIRE(channel.data.size() == 0);
}

TEST_CASE("Channel::add")
{
	Channel channel;
	channel.add(0.0);
	REQUIRE(channel.data.size() == 1);

	// clears data and only reserves storage of epos::NbrSAMP elems
	channel.init(epos::NbrSAMP);
	REQUIRE(channel.data.size() == 0);
}

TEST_CASE("Channel::getStats 1")
{
	Channel channel;
	Channel::Stats& stats = channel.stats;

	channel.add(0.0);
	channel.computeStats();
	REQUIRE(stats.min == 0.0);
	REQUIRE(stats.max == 0.0);
	REQUIRE(stats.idxmin == 0);
	REQUIRE(stats.idxmax == 0);
	REQUIRE(are_equal(stats.DC, 0.001));

	channel.add(1.0);
	channel.computeStats();
	REQUIRE(stats.min == 0.0);
	REQUIRE(stats.max == 1.0);
	REQUIRE(stats.idxmin == 0);
	REQUIRE(stats.idxmax == 1);

	channel.add(-1.0);
	channel.computeStats();
	REQUIRE(stats.min == -1.0);
	REQUIRE(stats.max == 1.0);
	REQUIRE(stats.idxmin == 2);
	REQUIRE(stats.idxmax == 1);
}

TEST_CASE("Channel::getStats 2")
{
	Channel ch;
	// clears data and only reserves storage of epos::NbrSAMP elems
	ch.init(epos::NbrSAMP);
	
	//Fill with 0
	for (int i = 0; i < epos::NbrSAMP; ++i)
		ch.add(0.0);

	REQUIRE(ch.data.size() == epos::NbrSAMP);
	REQUIRE(ch.data[0] == 0.0);

	ch.computeStats();
	REQUIRE(isInRange(ch.stats.DC, ch.stats.DC - ch.stats.DCstddev, ch.stats.DC + ch.stats.DCstddev));
	REQUIRE(ch.stats.min == 0.0);
	REQUIRE(ch.stats.max == 0.0);
	REQUIRE(ch.stats.idxmin == 0);
	REQUIRE(ch.stats.idxmax == 0);
}


TEST_CASE("Channel::getMean")
{
	Channel channel;
	Float mean;
	mean = channel.getMean(0);
	REQUIRE(mean == 0.0);
	mean = channel.getMean(10);
	REQUIRE(mean == 0.0);

	channel.add(2.0);
	mean = channel.getMean(1);
	REQUIRE(mean == 2.0);
	mean = channel.getMean(10);
	REQUIRE(mean == 2.0);

	channel.add(2.0);
	mean = channel.getMean(1);
	REQUIRE(mean == 2.0);
	mean = channel.getMean(2);
	REQUIRE(mean == 2.0);
	mean = channel.getMean(10);
	REQUIRE(mean == 2.0);
}

TEST_CASE("Channel::time")
{
	REQUIRE(Channel::time.size() == epos::NbrSAMP);
	for (int i = 0; i < epos::NbrSAMP; ++i) {
		REQUIRE(Channel::time[i] == Channel::getT(i));
	}
}
