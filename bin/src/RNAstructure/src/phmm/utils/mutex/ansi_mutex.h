#ifndef _ANSI_MUTEX_
#define _ANSI_MUTEX_

/*
Portable mutex, uses windows and pthread mutexes.
Refernces:
http://msdn.microsoft.com/en-us/library/ms686927%28v=VS.85%29.aspx
http://www.yolinux.com/TUTORIALS/LinuxTutorialPosixThreads.html#SYNCHRONIZATION
*/

#ifdef __unix__
// For linux multithreading
#include <pthread.h>
#include <errno.h>
#elif defined __APPLE__
#include <pthread.h>
#include <errno.h>
#elif defined _WIN32
#include <windows.h>
#else 
#error "The environment could not be determined, must specify manually."
#endif

class t_ansi_mutex
{
public:
	t_ansi_mutex();
	~t_ansi_mutex();

	// Declare the mutex variable, platform dependent, name is mutex.
#ifdef __unix__
	pthread_mutex_t mutex;
#elif defined __APPLE__
	pthread_mutex_t mutex;
#elif defined _WIN32
	HANDLE mutex;
#endif

	void lock_mutex();
	bool try_lock_mutex(); // Try locking the mutex.
	void release_mutex();

	// Free the resources related to the mutex.
	void free_mutex();
};

#endif // _ANSI_MUTEX_
