#include <stdio.h>
#include <stdlib.h>
#include "ansi_thread.h"

t_ansi_thread::t_ansi_thread(void* (*_thread_function)(void*), void* _thread_params)
{
	this->thread_function = _thread_function;
	this->thread_handle = NULL;
	this->id = ++thread_count;

	// Allocate and initialize thread parameters.
	this->thread_params = _thread_params;
	this->thread_state = THREAD_INITING;
}

t_ansi_thread::~t_ansi_thread()
{
}

int t_ansi_thread::thread_count = 0;

bool t_ansi_thread::run_thread()
{
#ifdef _WIN32
	// Create a Win32 thread.	
	DWORD dwGenericThread;
	char lszThreadParam[3];

	strcpy(lszThreadParam,"3");

	// Create the thread.
	this->thread_handle = CreateThread(NULL,
										0,
										(LPTHREAD_START_ROUTINE)this->thread_function,
										this->thread_params,
										0,
										&dwGenericThread);

	if(this->thread_handle == NULL)
	{
		// Just print the message that the client serving thread could not be created and return.
		DWORD last_error = GetLastError();
		printf("Could not create the socket serving thread @ %s(%d), last error is %d\n", __FILE__, __LINE__, (int)last_error);
		return(false);
	}
	else
	{
		return(true);
	}
#endif

#ifdef __unix__
	// Create pthread.
    pthread_attr_t attr;
 //   pthread_attr_init(&attr);
 //   pthread_attr_setstacksize (&attr, SOCKET_PTHREAD_STACK_SIZE);
	
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	int iret = pthread_create(&this->thread_handle, &attr, (void* (*)(void*))this->thread_function, (void*)this->thread_params);

	if(iret != 0)
	{
		printf("Could not create the socket serving thread @ %s(%d)\n", __FILE__, __LINE__);
		return(false);
	}
	else
	{
		// Detach the thread from the process. (For resource freeing)
		return(true);
	}
#endif

#ifdef __APPLE__
        // Create pthread.
	pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        int iret = pthread_create(&this->thread_handle, &attr, (void* (*)(void*))this->thread_function, (void*)this->thread_params);

        if(iret != 0)
	  {
	    printf("Could not create the socket serving thread @ %s(%d)\n", __FILE__, __LINE__);
	    return(false);
	  }
        else
	  {
	    // Detach the thread from the process. (For resource freeing)
	    return(true);
	  }
#endif
}

void t_ansi_thread::wait_thread()
{
#ifdef _WIN32
	WaitForSingleObject(this->thread_handle, INFINITE);
#endif 

#ifdef __unix__
	// Wait until thread finishes.
	//while(this->thread_state != THREAD_TERMINAL){};
	void* status = NULL;
	int rc = pthread_join(this->thread_handle, &status);
	if (rc) 
	{
		printf("ERROR; return code from pthread_join() is %d\n", rc);
		exit(-1);
	}
#endif

#ifdef __APPLE__
        // Wait until thread finishes.
	void* status = NULL;
        int rc = pthread_join(this->thread_handle, &status);
        if (rc)
	  {
	    printf("ERROR; return code from pthread_join() is %d\n", rc);
	    exit(-1);
	  }
#endif
}

THREAD_HANDLE t_ansi_thread::current_thread_handle() {
#ifdef _WIN32
	return GetCurrentThread(); // see also GetCurrentThreadId();
#endif 
#ifdef __unix__
    return pthread_self();
#endif
#ifdef __APPLE__
   return pthread_self();
#endif
}

void t_ansi_thread::exit_thread()
{
	this->thread_state = THREAD_TERMINAL;

#ifdef __unix__
        pthread_exit(NULL);
#endif

#ifdef __APPLE__
        pthread_exit(NULL);
#endif
}

