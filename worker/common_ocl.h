#ifndef COMMON_OCL
#define COMMON_OCL

#define _WORK_ITEMS_PER_DIMENSION_  16

#include <stdio.h>
#include <string.h>
#include <CL/cl.h>


/* Setting opencl by inititalizing 1 device on 1 context with num_queues queues */
void set_opencl_env_multiple_queues(int num_queues, cl_command_queue **list_queues, cl_context *context);

/* Releases the opencl queues and memory used */
void release_opencl(int num_queues, cl_command_queue **list_queues, cl_context *context);

/**
 * Builds a kernel for a given device and context from a source file.-
 *
 */
void build_kernel_from_file (cl_context *context, 
                             char *kernel_name, 
                             const char *file_name, 
                             cl_kernel *kernel);

/**
 * Builds a kernel for a given device and context from source. 
 *
 */
void build_kernel (cl_context *context, 
                   char *kernel_name, 
                   const char **kernel_source, 
                   cl_kernel *kernel);

/**
 * Attaches the received parameter as a kernel argument.-
 *
 */
void set_kernel_arg (cl_kernel *kernel, 
                     cl_uint arg_index, 
                     size_t arg_size, 
                     const void *arg_value);

/**
 * Marks a kernel parameter as local memory.-
 *
 */
void set_local_mem (cl_kernel *kernel, cl_uint arg_index, size_t arg_size);

/**
 * Creates an OpenCL buffer, returning a pointer to it in the last parameter.-
 *
 */
void create_buffer (cl_context *context, 
                    cl_mem_flags flags, 
                    size_t size, 
                    void *host_ptr, 
                    cl_mem *dev_ptr);

/**
 * Reads an OpenCL buffer from the device.-
 *
 */
void read_buffer (cl_command_queue *queue, 
                  cl_mem *dev_ptr, 
                  cl_bool blocking_read, 
                  size_t size, 
                  void *host_ptr);

/**
 * Writes an OpenCL buffer to the device.-
 *
 */
void write_buffer (cl_command_queue *queue, 
                   cl_mem *dev_ptr, 
                   cl_bool blocking_write, 
                   size_t size, 
                   const void *host_ptr);

/**
 * Executes the specified kernel.-
 *
 */
void run_kernel (cl_command_queue *queue, cl_kernel *kernel);

/**
 * Check and prints OpenCL errors.-
 *
 */
void check_error (int status, char *msg);

/* Internal helper functions not meant to be directly used within an application */
void print_device_information(cl_device_id *device);
void get_amd_platform(cl_platform_id *platform);
int  set_device_and_context(cl_platform_id *platform, cl_device_id *device, cl_context *context);

#endif /* COMMON_OCL */
