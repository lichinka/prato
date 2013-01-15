#include "ocl_common.h"



/**
 * Reads the kernel source file from and saves it in the output parameter,
 * returning the number of bytes read.
 *
 * file_name        The file from which to read the kernel source.
 * kernel_source    output parameter, where the file contents are saved.-
 *
 */
static int read_kernel_file (const char *file_name, 
                             const char *kernel_source)
{
    char *read_ptr = (char *) kernel_source;
    FILE *fr = fopen (file_name, "r");

    if (fr == NULL)
    {
        fprintf (stderr, "ERROR Cannot open [%s] for reading\n", file_name);
        exit (1);
    }

    while (fgets (read_ptr, 1024, fr) != NULL)
       read_ptr += strlen (read_ptr);

    fclose (fr);

    return (strlen (kernel_source) + 1);
}



/**
 * Displays a friendly representation of the OpenCl error codes 
 * to the standard error stream.-
 *
 */
static void print_opencl_error (int error)
{
    switch (error)
    {
        case (CL_BUILD_PROGRAM_FAILURE):
            fprintf (stderr, "Build program failure\n");
            break;
        case (CL_COMPILER_NOT_AVAILABLE):
            fprintf (stderr, "Compiler not available\n");
            break;
        case (CL_DEVICE_NOT_AVAILABLE):
            fprintf (stderr, "Device not available\n");
            break;
        case (CL_INVALID_ARG_SIZE):
            fprintf (stderr, "If the argument is not a memory object, [arg_size] does not match the size of the data type.\n");
            fprintf (stderr, "If the argument is a memory object, [arg_size] != sizeof(cl_mem).\n");
            fprintf (stderr, "If the argument is declared with the __local qualifier, [arg_size] is 0.\n");
            break;
        case (CL_INVALID_CONTEXT):
            fprintf (stderr, "Invalid context\n");
            break;
        case (CL_INVALID_BINARY):
            fprintf (stderr, "Invalid binary\n");
            break;
        case (CL_INVALID_BUFFER_SIZE):
            fprintf (stderr, "Invalid buffer size\n");
            break;
        case (CL_INVALID_BUILD_OPTIONS):
            fprintf (stderr, "Invalid build options\n");
            break;
        case (CL_INVALID_DEVICE):
            fprintf (stderr, "Invalid device\n");
            break;
        case (CL_INVALID_EVENT):
            fprintf (stderr, "Invalid event\n");
            break;
        case (CL_INVALID_HOST_PTR):
            fprintf (stderr, "Invalid host pointer\n");
            break;
        case (CL_INVALID_KERNEL_ARGS):
            fprintf (stderr, "Kernel argument values have not been specified\n");
            break;
        case (CL_INVALID_MEM_OBJECT):
            fprintf (stderr, "Memory objects are not valid or are not buffer objects\n");
            break;
        case (CL_INVALID_OPERATION):
            fprintf (stderr, "Invalid operation\n");
            break;
        case (CL_INVALID_PLATFORM):
            fprintf (stderr, "Invalid platform\n");
            break;
        case (CL_INVALID_PROGRAM):
            fprintf (stderr, "Invalid program\n");
            break;
        case (CL_INVALID_PROPERTY):
            fprintf (stderr, "Invalid property\n");
            break;
        case (CL_INVALID_VALUE):
            fprintf (stderr, "Invalid value\n");
            break;
        case (CL_MEM_OBJECT_ALLOCATION_FAILURE):
            fprintf (stderr, "Cl_mem object allocation failure\n");
            break;
        case (CL_OUT_OF_HOST_MEMORY):
            fprintf (stderr, "Out of host memory\n");
            break;
        case (CL_OUT_OF_RESOURCES):
            fprintf (stderr, "Out of resources\n");
            break;
        case (CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST):
            fprintf (stderr, "invalid context\n");
            break;
        default:
            fprintf (stderr, "Unkown error code [%d]\n", error);
    }
}



/**
 * Checks that the structure passed has been correctly initialized.-
 */
static void check_is_initialized (OCL_objects *ocl_obj)
{
    if (ocl_obj->is_initialized == 0)
    {
        fprintf (stderr,
                 "ERROR OpenCL has not been correctly initialized\n");
        exit (1);
    }
}



/**
 * Checks that the command queue with this index exists.-
 */
static void check_queue (OCL_objects *ocl_obj,
                         int queue_index)
{
    if (queue_index >= ocl_obj->queue_count)
    {
        fprintf (stderr,
                 "ERROR Command queue with index (%d) does not exist\n",
                 queue_index);
        exit (1);
    }
    /*
    else
    {
        cl_device_id *param_value = (cl_device_id *) malloc (sizeof (cl_device_id));
        clGetCommandQueueInfo (ocl_obj->queues[queue_index],
                               CL_QUEUE_DEVICE,
                               sizeof (cl_device_id),
                               param_value,
                               NULL);
        printf ("INFO command queue\n");
        print_device_information (param_value);
        free (param_value);
    }
    */
}



void check_error (int status, char *msg)
{
    if (status != CL_SUCCESS)
    {
        fprintf (stderr, "*** OpenCL ERROR at %s: ", msg);
        print_opencl_error (status);
    }
}



void print_device_information (cl_device_id *device)
{
  cl_device_type  device_type;
  char            device_name_buffer[100];  
  char            string_gpu[] = "GPU";
  char            string_cpu[] = "CPU";
  char           *string_device_type;

  clGetDeviceInfo(*device, CL_DEVICE_TYPE, sizeof(device_type), &device_type, NULL);
  string_device_type = string_gpu;
  if(device_type == CL_DEVICE_TYPE_CPU) { string_device_type = string_cpu; }

  clGetDeviceInfo(*device, CL_DEVICE_NAME, sizeof(device_name_buffer), &device_name_buffer, NULL);
  printf("%s device: %s\n", string_device_type, device_name_buffer);
}



/**
 * Initializes the OpenCL platform. 
 * This function should be called before any other one.
 *
 * ocl_obj          A pointer to the uninitialized OCL_objects structure on
 *                  the user's side;
 * number_of_queues the number of queues to initialize.-
 *
 */
void init_opencl (OCL_objects *ocl_obj, int number_of_queues)
{
    int i;

    if (ocl_obj->is_initialized == 1)
        fprintf (stderr,
                 "WARNING OpenCL has already been initialized [%d]\n",
                 ocl_obj->is_initialized);

    get_platform (&(ocl_obj->platform_id));
    set_device_and_context (&(ocl_obj->platform_id),
                            &(ocl_obj->device_id),
                            &(ocl_obj->context));
    ocl_obj->queues = (cl_command_queue *) calloc (number_of_queues,
                                                   sizeof (cl_command_queue));
    ocl_obj->events = (cl_event *) calloc (number_of_queues,
                                           sizeof (cl_event));
    for (i = 0; i < number_of_queues; i ++)
    {
        ocl_obj->queues[i] = clCreateCommandQueue (ocl_obj->context, 
                                                   ocl_obj->device_id, 
                                                   0,
                                                   &(ocl_obj->status));
        check_error (ocl_obj->status,
                     "Initialize OpenCL: create command queue");
    }
    ocl_obj->queue_count = number_of_queues;
    ocl_obj->is_initialized = 1;
}



/**
 * Initializes the OpenCL platform. 
 * This function should be called before any other one.-
 *
 */
void init_platform (cl_platform_id *platform, 
                    cl_device_id *device, 
                    cl_context *context,
                    cl_command_queue *queue)
{
    cl_int status;
    get_platform (platform);
    set_device_and_context (platform,
                            device,
                            context);
    *queue = clCreateCommandQueue (*context, 
                                   *device, 
                                   0, 
                                   &status);
    check_error (status, "Initialize platform");
}



/**
 * Queries the system for a OpenCL platform, returning it
 * in the input parameter.-
 */
void get_platform (cl_platform_id *platform)
{
	int  i;
	int  status;
	char vendor_buffer [512];

	cl_platform_id *list_platforms;

	uint num_platforms;
	status = clGetPlatformIDs(0, NULL, &num_platforms);
	check_error (status, "Get platform ID");

	if (num_platforms > 0)
	{
		list_platforms = malloc (num_platforms * sizeof (cl_platform_id));
		status = clGetPlatformIDs (num_platforms, list_platforms, NULL);
		check_error (status, "Get platform ID");

		for (i = 0; i < num_platforms; i++)
		{
			status = clGetPlatformInfo (list_platforms[i], 
				  						CL_PLATFORM_VENDOR, 
				  						sizeof (vendor_buffer), 
				  						vendor_buffer, 
				  						NULL);
			check_error (status, "Get platform info");

			//
			// save the first valid platform
			//
			*platform = list_platforms[i];
			break;
		}
		free (list_platforms);
	}
	else
		fprintf (stderr, "*** ERROR: No OpenCL platforms found\n");
}



int set_device_and_context (cl_platform_id *platform, 
                            cl_device_id *device, 
                            cl_context *context)
{
  int status;

  uint          num_devices;
  cl_device_id *list_devices;

  uint           device_number = 0; /* Get first device it sees */
  cl_device_type device_type   = CL_DEVICE_TYPE_GPU;

  status = clGetDeviceIDs (*platform, device_type, 0, NULL, &num_devices);
  if (status == CL_DEVICE_NOT_FOUND)
  {
      fprintf (stderr,
               "GPU not found. Falling back to CPU.\n");
      device_type = CL_DEVICE_TYPE_CPU;
  }
  check_error (status, "Get device ID");

  if (num_devices > 0)
  {
    list_devices = malloc(num_devices * sizeof(cl_device_id));

    /* Set a device */
    status = clGetDeviceIDs(*platform, device_type, num_devices, list_devices, NULL);
    check_error(status, "Get device ID");
    *device = list_devices[device_number];

    print_device_information(device);

    free(list_devices);
  } else {
    printf("No device available.\n");
    return 0;
  }

  /* Set a context */
  *context = clCreateContext(NULL, 1, device, NULL, NULL, &status);
  check_error(status, "Create Context");

  return 1;
}



void get_device_from_context(cl_context *context, cl_device_id *device)
{
  uint          num_devices;
  cl_device_id *list_devices;

  clGetContextInfo(*context, CL_CONTEXT_NUM_DEVICES, sizeof(uint), &num_devices, NULL);

  list_devices = malloc(num_devices * sizeof(cl_device_id));

  clGetContextInfo(*context, CL_CONTEXT_DEVICES, num_devices * sizeof(cl_device_id), list_devices, NULL);

  *device = list_devices[0];

  free(list_devices);
}



/**
 * Sets N queues for 1 device using 1 context from 1 platform.
 */
void set_opencl_env_multiple_queues(int num_queues, cl_command_queue **list_queues, cl_context *context)
{
  int i;
  int status;

  cl_platform_id platform;
  cl_device_id   device;

  get_platform (&platform);

  set_device_and_context(&platform, &device, context);

  *list_queues = malloc(num_queues * sizeof(cl_command_queue));

 
  for (i = 0; i < num_queues; i++)
  {
    (*list_queues)[i] = clCreateCommandQueue(*context, device, 0, &status);
    check_error(status, "Create command queue");
  }
}



/**
 * Attaches the received value as a kernel argument.-
 *
 * ocl_obj          A pointer to the initialized OCL_objects structure on the
 *                  user's side;
 * arg_index        argument index in the kernel-function prototype;
 * arg_size         size (in bytes) of the argument being passed;
 * arg_value_ptr    pointer to the argument value.-
 *
 */
void set_kernel_value_arg (OCL_objects *ocl_obj,
                           cl_uint arg_index, 
                           size_t arg_size, 
                           const void *arg_value_ptr)
{
    check_is_initialized (ocl_obj);
    char msg [1024];

    snprintf (msg, 
              1024, 
              "Set kernel value argument with index [%d]", 
              arg_index);
    ocl_obj->status = clSetKernelArg (ocl_obj->kernel,
                                      arg_index,
                                      arg_size,
                                      arg_value_ptr);
    check_error (ocl_obj->status, 
                 msg);
}



/**
 * Attaches the received memory object as a kernel argument.-
 *
 * ocl_obj          A pointer to the initialized OCL_objects structure on the
 *                  user's side;
 * arg_index        argument index in the kernel-function prototype;
 * arg_value_ptr    pointer to the argument value.-
 *
 */
void set_kernel_mem_arg (OCL_objects *ocl_obj,
                         cl_uint arg_index, 
                         const void *arg_value_ptr)
{
    check_is_initialized (ocl_obj);
    char msg [1024];

    snprintf (msg, 
              1024, 
              "Set kernel memory argument with index [%d]", 
              arg_index);
    ocl_obj->status = clSetKernelArg (ocl_obj->kernel,
                                      arg_index,
                                      sizeof (cl_mem),
                                      arg_value_ptr);
    check_error (ocl_obj->status, 
                 msg);
}



/**
 * Attaches the received integer value as a kernel argument.-
 *
 * ocl_obj          A pointer to the initialized OCL_objects structure on the
 *                  user's side;
 * arg_index        argument index in the kernel-function prototype;
 * arg_value_ptr    pointer to the argument value.-
 *
 */
void set_kernel_int_arg (OCL_objects *ocl_obj,
                         cl_uint arg_index, 
                         const int *arg_value_ptr)
{
    set_kernel_value_arg (ocl_obj,
                          arg_index,
                          sizeof (*arg_value_ptr),
                          arg_value_ptr);
}



/**
 * Attaches the received float value as a kernel argument.-
 *
 * ocl_obj          A pointer to the initialized OCL_objects structure on the
 *                  user's side;
 * arg_index        argument index in the kernel-function prototype;
 * arg_value_ptr    pointer to the argument value.-
 *
 */
void set_kernel_float_arg (OCL_objects *ocl_obj,
                           cl_uint arg_index, 
                           const float *arg_value_ptr)
{
    set_kernel_value_arg (ocl_obj,
                          arg_index,
                          sizeof (*arg_value_ptr),
                          arg_value_ptr);
}



/**
 * Attaches the received double value as a kernel argument.-
 *
 * ocl_obj          A pointer to the initialized OCL_objects structure on the
 *                  user's side;
 * arg_index        argument index in the kernel-function prototype;
 * arg_value_ptr    pointer to the argument value.-
 *
 */
void set_kernel_double_arg (OCL_objects *ocl_obj,
                           cl_uint arg_index, 
                           const double *arg_value_ptr)
{
    set_kernel_value_arg (ocl_obj,
                          arg_index,
                          sizeof (*arg_value_ptr),
                          arg_value_ptr);
}



/**
 * Marks a kernel parameter as local memory.-
 *
 * ocl_obj          A pointer to the initialized OCL_objects structure on the
 *                  user's side;
 * arg_index        argument index in the kernel-function prototype;
 * arg_size         size (in bytes) of the allocated local memory.-
 *
 */
void set_local_mem (OCL_objects *ocl_obj,
                    cl_uint arg_index, 
                    size_t arg_size)
{
    check_is_initialized (ocl_obj);
    char msg [1024];

    snprintf (msg, 
              1024, 
              "Set kernel argument with index [%d] as local memory, with size [%ld]",
              arg_index,
              arg_size);

    ocl_obj->status = clSetKernelArg (ocl_obj->kernel,
                                      arg_index,
                                      arg_size,
                                      NULL);
    check_error (ocl_obj->status, 
                 msg);
}



/**
 * Creates and returns a OpenCL buffer, i.e. cl_mem object.
 *
 * ocl_obj      A pointer to the initialized OCL_objects structure on the
 *              user's side;
 * flags        memory flags, as defined by OpenCL;
 * size         size (in bytes) of the buffer to create.-
 *
 */
cl_mem create_buffer (OCL_objects *ocl_obj,
                      cl_mem_flags flags, 
                      size_t size)
{
    check_is_initialized (ocl_obj);
    cl_mem ret_value = clCreateBuffer (ocl_obj->context,
                                       flags,
                                       size,
                                       NULL,
                                       &(ocl_obj->status));
    check_error (ocl_obj->status, 
                 "Create buffer");
    return ret_value;
}



/**
 * Reads an OpenCL buffer from the device, waiting for the operation to
 * finish before returning.
 *
 * ocl_obj      A pointer to the initialized OCL_objects structure on the
 *              user's side;
 * queue_index  index of the command queue on which the read is performed;
 * cl_mem_ptr   pointer to a valid cl_mem object, from which to copy;
 * size         size (in bytes) of the data to read from the device;
 * target_ptr   pointer to the data to which to copy, on the host.-
 *
 */
void read_buffer_blocking (OCL_objects *ocl_obj,
                           int queue_index,
                           cl_mem *cl_mem_ptr, 
                           size_t size, 
                           void *target_ptr)
{
    check_is_initialized (ocl_obj);
    check_queue (ocl_obj, 
                 queue_index);
    ocl_obj->status = clEnqueueReadBuffer (ocl_obj->queues[queue_index],
                                           *cl_mem_ptr,
                                           CL_TRUE,
                                           0,
                                           size,
                                           target_ptr,
                                           0,
                                           NULL,
                                           &(ocl_obj->events[queue_index]));
    check_error (ocl_obj->status, 
                 "Read buffer");
}



/**
 * Enqueues a write operation of an OpenCL buffer to the device, returning
 * the associated event object.
 *
 * ocl_obj      A pointer to the initialized OCL_objects structure on the
 *              user's side;
 * queue_index  index of the command queue on which the write is performed;
 * cl_mem_ptr   pointer to a valid cl_mem object;
 * size         size (in bytes) of the data to write to the device;
 * source_ptr   to the data from which to copy, on the host.-
 *
 */
cl_event* write_buffer (OCL_objects *ocl_obj,
                        int queue_index,
                        cl_mem *cl_mem_ptr, 
                        size_t size, 
                        const void *source_ptr)
{
    check_is_initialized (ocl_obj);
    check_queue (ocl_obj, 
                 queue_index);
    ocl_obj->status = clEnqueueWriteBuffer (ocl_obj->queues[queue_index],
                                            *cl_mem_ptr,
                                            CL_TRUE,
                                            0,
                                            size,
                                            source_ptr,
                                            0,
                                            NULL,
                                            &(ocl_obj->events[queue_index]));
    check_error (ocl_obj->status, 
                 "Write buffer - non blocking");
    return &(ocl_obj->events[queue_index]);
}



/**
 * Writes an OpenCL buffer to the device, waiting for the operation to
 * finish before returning.
 *
 * ocl_obj      A pointer to the initialized OCL_objects structure on the
 *              user's side;
 * queue_index  index of the command queue on which the write is performed;
 * cl_mem_ptr   pointer to a valid cl_mem object;
 * size         size (in bytes) of the data to write to the device;
 * source_ptr   to the data from which to copy, on the host.-
 *
 */
void write_buffer_blocking (OCL_objects *ocl_obj,
                            int queue_index,
                            cl_mem *cl_mem_ptr, 
                            size_t size, 
                            const void *source_ptr)
{
    check_is_initialized (ocl_obj);
    check_queue (ocl_obj, 
                 queue_index);
    ocl_obj->status = clEnqueueWriteBuffer (ocl_obj->queues[queue_index],
                                            *cl_mem_ptr,
                                            CL_TRUE,
                                            0,
                                            size,
                                            source_ptr,
                                            0,
                                            NULL,
                                            NULL);
    check_error (ocl_obj->status, 
                 "Write buffer - blocking");
}



/**
 * Executes the kernel, using the received 2D range, waiting for it to finish.
 *
 * ocl_obj          A pointer to the initialized OCL_objects structure on the
 *                  user's side;
 * queue_index      index of the command queue on which the kernel execution is 
 *                  performed;
 * global_offsets   offsets used for the kernel range (2 elements). If NULL
 *                  is given, the offsets used are zero;
 * global_sizes     global execution range sizes (2 elements);
 * local_sizes      local execution range sizes (2 elements);
 *
 */
void run_kernel_2D_blocking (OCL_objects *ocl_obj,
                             int queue_index,
                             const size_t *global_offsets,
                             const size_t *global_sizes,
                             const size_t *local_sizes)
{
    size_t goff [] = {0, 0};
    check_is_initialized (ocl_obj);
    check_queue (ocl_obj,
                 queue_index);
    if (global_offsets != NULL)
    {
        goff[0] = global_offsets[0];
        goff[1] = global_offsets[1];
    }
    ocl_obj->status = clEnqueueNDRangeKernel (ocl_obj->queues[queue_index],
                                              ocl_obj->kernel,
                                              2,
                                              goff,
                                              global_sizes,
                                              local_sizes,
                                              0,
                                              NULL,
                                              &(ocl_obj->events[queue_index]));
    check_error (ocl_obj->status, 
                 "Run kernel with 2D range");
    clWaitForEvents (1,
                     &(ocl_obj->events[queue_index]));
}



/**
 * Releases OpenCL, deallocating all referenced memory.
 */
void release_opencl (int num_queues, 
                     cl_command_queue **list_queues, 
                     cl_context *context)
{
  int i;

  // Flushes and finishes all queues in the list of queues 
  for (i = 0; i < num_queues; i++)
  {
    clFinish((*list_queues)[i]);
  }

  clReleaseContext(*context);

  // Free memory used by list of pointers
  free(*list_queues);
}



/**
 * Deactivates OpenCL, deallocating all referenced memory.
 *
 * ocl_obj      A pointer to the initialized OCL_objects structure on the
 *              user's side.-
 *
 */
void deactivate_opencl (OCL_objects *ocl_obj)
{
    check_is_initialized (ocl_obj);
    release_opencl (ocl_obj->queue_count,
                    &(ocl_obj->queues),
                    &(ocl_obj->context));

    // Free memory used by the list of events
    free (ocl_obj->events);
}



/**
 * Builds a kernel for a given device and context from a source file.-
 *
 * ocl_obj      A pointer to the initialized OCL_objects structure on the
 *              user's side;
 * kernel_name  name of the kernel to activate;
 * file_name    file containing the OpenCL kernel code;
 * options      compile-time options, like include directories or defines.-
 *
 */
void build_kernel_from_file (OCL_objects *ocl_obj,
                             char *kernel_name, 
                             const char *file_name, 
                             const char *options)
{
    check_is_initialized (ocl_obj);
    int bytes_read;
    char *source = (char *) calloc (1024000,
                                    sizeof (char));

    bytes_read = read_kernel_file (file_name, 
                                   source);
    if (bytes_read < 1)
    {
        fprintf (stderr, 
                 "ERROR Could not load kernel source from <%s>",
                 file_name);
        exit (1);
    }
    build_kernel (&(ocl_obj->context),
                  kernel_name,
                  (const char **) &source,
                  options,
                  &(ocl_obj->kernel));
    free (source);
}



/**
 * Builds a kernel for a given device and context from string buffer.-
 *
 */
void build_kernel (cl_context *context, 
                   char *kernel_name, 
                   const char **kernel_source, 
                   const char *options,
                   cl_kernel *kernel)
{
    int          status;
    cl_device_id device;
    uint         num_devices = 1;
    uint         num_strings_in_kernel_source = 1;

    get_device_from_context (context, &device);
    cl_program program = clCreateProgramWithSource (*context, 
                                                    num_strings_in_kernel_source, 
                                                    kernel_source, 
                                                    NULL, 
                                                    &status);
    check_error (status, "Create OCL program with source");
    int build_status = clBuildProgram (program, 
                                       num_devices, 
                                       &device, 
                                       options,
                                       NULL, 
                                       NULL);
    check_error (status, "Build OCL program");
    if (build_status != CL_SUCCESS)
    {
        int error_buffer_size = 102400;
        char error_buffer [error_buffer_size];

        clGetProgramBuildInfo (program, 
                               device, 
                               CL_PROGRAM_BUILD_LOG, 
                               error_buffer_size, 
                               error_buffer, 
                               NULL);
        fprintf (stderr, 
                 "*** Compilation failed. Log follows:\n\n%s\n", 
                 error_buffer);
        exit (1);
    }
    *kernel = clCreateKernel (program, 
                              kernel_name, 
                              NULL);
}


/* EOF */
