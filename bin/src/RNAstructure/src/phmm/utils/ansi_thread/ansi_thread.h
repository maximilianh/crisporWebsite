#ifndef _ANSI_THREAD_
#define _ANSI_THREAD_

/*
t_socket_thread class: Represents a thread for any type of connection oriented socket, this can be on either side of connection, i.e., it can be
used in client or server. 
*/

enum{THREAD_INITING, THREAD_RUNNING, THREAD_TERMINAL}; 

#ifdef __unix__
// For linux multithreading
	#include <pthread.h>
	#define THREAD_HANDLE pthread_t
#elif defined __APPLE__
	#include <pthread.h>
	#define THREAD_HANDLE pthread_t
#elif defined _WIN32
	#include <windows.h>
	#define THREAD_HANDLE HANDLE
#else 
	#error "The environment could not be determined, must specify manually."
#endif

class t_ansi_thread
{
public:
	static int thread_count;
	t_ansi_thread(void* (*_thread_function)(void*), void* _thread_params);
	~t_ansi_thread();

	void* (*thread_function)(void*);

	void* thread_params;

	// State of the thread, necessary to keep track of detached pthread's. 
	volatile int thread_state; 

	int id;

	bool run_thread(); // Start the socket thread. (Start protocol_implementing_function as a thread)
	void wait_thread(); // Wait for the socket thread to finish.
	void exit_thread(); // exit this thread.
	//! Get the ID/HANDLE of the this thread object.
	inline THREAD_HANDLE get_handle() { return thread_handle; }
	//! Get the ID/HANDLE of the currently executing thread (without any reference to a thread object)
	static THREAD_HANDLE current_thread_handle();

	private:
	//! In Windows, this is a Real (Win32) handle of accepted client.
	//!    Use WaitForSingleObject to wait for the threads to finish.
	//! In Linux and OSX, it is a pthread handle of accepted client.
	//!    Use pthread_join (http://www.yolinux.com/TUTORIALS/LinuxTutorialPosixThreads.html#BASICS)
	//!    to wait for the threads to finish.
	THREAD_HANDLE thread_handle;
};

#endif // _ANSI_THREAD_
