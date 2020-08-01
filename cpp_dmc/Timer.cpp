
#include "Timer.h"

cpp_mc::Timer::Timer()
{
	start = std::chrono::high_resolution_clock::now();
	stop = std::chrono::high_resolution_clock::now();
}

double cpp_mc::Timer::GetMilisecondsElapsed()
{
	if (isrunning)
	{
		auto elapsed = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start);
		return elapsed.count();
	}
	else
	{
		auto elapsed = std::chrono::duration<double, std::milli>(stop - start);
		return elapsed.count();
	}
}

void cpp_mc::Timer::Restart()
{
	isrunning = true;
	start = std::chrono::high_resolution_clock::now();
}

bool cpp_mc::Timer::Stop()
{
	if (!isrunning)
		return false;
	else
	{
		stop = std::chrono::high_resolution_clock::now();
		isrunning = false;
		return true;
	}
}

bool cpp_mc::Timer::Start()
{
	if (isrunning)
	{
		return false;
	}
	else
	{
		start = std::chrono::high_resolution_clock::now();
		isrunning = true;
		return true;
	}
}
