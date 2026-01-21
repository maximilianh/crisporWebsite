#include <stdio.h>
#include <stdlib.h>
#include "ansi_mutex.h"

t_ansi_mutex::t_ansi_mutex()
{
#ifdef __unix__
	pthread_mutex_init(&this->mutex, NULL);
#elif defined __APPLE__
	pthread_mutex_init(&this->mutex, NULL);
#elif defined _WIN32
	this->mutex = CreateMutex(NULL,              // default security attributes
								FALSE,             // initially not owned
								NULL);             // unnamed mutex

    if (this->mutex == NULL) 
    {
        printf("CreateMutex error: %d\n", GetLastError());
        exit(0);
    }
#endif
}

t_ansi_mutex::~t_ansi_mutex()
{
	this->free_mutex();
}

void t_ansi_mutex::free_mutex()
{
#ifdef __unix__
	// Do nothing?
#elif defined __APPLE__
  // Do nothing?
#elif defined _WIN32
	CloseHandle(this->mutex);
#endif
}

// Following is blocking: Blocks until the mutex is released.
void t_ansi_mutex::lock_mutex()
{
#ifdef __unix__
	pthread_mutex_lock(&this->mutex);
#elif defined __APPLE__
	pthread_mutex_lock(&this->mutex);
#elif defined _WIN32
	// WaitForSingleObject wait for the mutex and gets the ownership when mutex is released.
	DWORD dwWaitResult = WaitForSingleObject( 
				this->mutex,    // handle to mutex
				INFINITE);  // no time-out interval
#endif
}

bool t_ansi_mutex::try_lock_mutex()
{
#ifdef __unix__
	int ret = pthread_mutex_trylock(&this->mutex);
	if(ret == 0)
	{
		return(true);
	}
	else if(ret == EBUSY)
	{
		return(false);
	}
	else if(ret == EINVAL )
	{
		printf("try_lock failed @ %s(%d): Uninitialized mutex.\n", __FILE__, __LINE__);
		exit(0);
	}
	else if(ret == EFAULT)
	{
		printf("try_lock failed @ %s(%d): Invalid mutex pointer.\n", __FILE__, __LINE__);
		exit(0);
	}
	else
	{
		printf("try_lock failed @ %s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}
#elif defined __APPLE__
        int ret = pthread_mutex_trylock(&this->mutex);
        if(ret == 0)
	  {
	    return(true);
	  }
        else if(ret == EBUSY)
	  {
	    return(false);
	  }
        else if(ret == EINVAL )
	  {
	    printf("try_lock failed @ %s(%d): Uninitialized mutex.\n", __FILE__, __LINE__);
	    exit(0);
	  }
        else if(ret == EFAULT)
	  {
	    printf("try_lock failed @ %s(%d): Invalid mutex pointer.\n", __FILE__, __LINE__);
	    exit(0);
	  }
        else
	  {
	    printf("try_lock failed @ %s(%d)\n", __FILE__, __LINE__);
	    exit(0);
	  }
#elif defined _WIN32
	// WaitForSingleObject wait for the mutex and gets the ownership when mutex is released.
	DWORD dwWaitResult = WaitForSingleObject( 
				this->mutex,    // handle to mutex
				0);  // no time-out interval

	if(dwWaitResult == WAIT_OBJECT_0)
	{
		return(true);
	}
	else if(dwWaitResult == WAIT_ABANDONED)
	{
		return(true);
	}
	else if(dwWaitResult == WAIT_TIMEOUT)
	{
		return(false);
	}
	else if(dwWaitResult == WAIT_FAILED)
	{
		printf("try_lock failed @ %s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}
#endif
	return true;
}


void t_ansi_mutex::release_mutex()
{
#ifdef __unix__
	pthread_mutex_unlock(&this->mutex);
#elif defined __APPLE__
	pthread_mutex_unlock(&this->mutex);
#elif defined _WIN32
    // Release ownership of the mutex object
    if (!ReleaseMutex(this->mutex)) 
    { 
        printf("Could not release mutex: %d\n", GetLastError());
        exit(0);
    } 
#endif
}

