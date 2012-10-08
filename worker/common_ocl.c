#include "worker/common_ocl.h"




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
    if (error == CL_BUILD_PROGRAM_FAILURE)  fprintf (stderr, "Build program failure");
    if (error == CL_COMPILER_NOT_AVAILABLE) fprintf (stderr, "Compiler not available");
    if (error == CL_DEVICE_NOT_AVAILABLE)   fprintf (stderr, "Device not available");
    if (error == CL_INVALID_CONTEXT)        fprintf (stderr, "Invalid context");
    if (error == CL_INVALID_BINARY)         fprintf (stderr, "Invalid binary");
    if (error == CL_INVALID_BUFFER_SIZE)    fprintf (stderr, "Invalid buffer size");
    if (error == CL_INVALID_BUILD_OPTIONS)  fprintf (stderr, "Invalid build options");
    if (error == CL_INVALID_DEVICE)         fprintf (stderr, "Invalid device");
    if (error == CL_INVALID_EVENT)          fprintf (stderr, "Invalid event");
    if (error == CL_INVALID_HOST_PTR)       fprintf (stderr, "Invalid host pointer");
    if (error == CL_INVALID_OPERATION)      fprintf (stderr, "Invalid operation");
    if (error == CL_INVALID_PLATFORM)       fprintf (stderr, "Invalid platform");
    if (error == CL_INVALID_PROGRAM)        fprintf (stderr, "Invalid program");
    if (error == CL_INVALID_PROPERTY)       fprintf (stderr, "Invalid property");
    if (error == CL_INVALID_VALUE)          fprintf (stderr, "Invalid value");
    if (error == CL_MEM_OBJECT_ALLOCATION_FAILURE) fprintf (stderr, "Cl_mem object allocation failure");
    if (error == CL_OUT_OF_HOST_MEMORY)     fprintf (stderr, "Out of host memory");
    if (error == CL_OUT_OF_RESOURCES)       fprintf (stderr, "Out of resources");
    if (error == CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST) fprintf (stderr, "invalid context");
    fprintf (stderr, "\n");
}



void check_error (int status, char *msg)
{
  if (status != CL_SUCCESS)
  {
    fprintf (stderr, "OpenCL error at %s: ", msg);
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



void get_amd_platform(cl_platform_id *platform)
{
  int  i;
  int  status;
  char vendor_buffer[100];
  char amd_vendor[] = "Advanced Micro Devices, Inc.";  

  cl_platform_id *list_platforms;

  uint num_platforms;
  status = clGetPlatformIDs(0, NULL, &num_platforms);
  check_error(status, "Get platform ID");

  if (num_platforms > 0)
  {
    list_platforms = malloc(num_platforms * sizeof(cl_platform_id));
    status = clGetPlatformIDs(num_platforms, list_platforms, NULL);
    check_error(status, "Get platform ID");

    for (i = 0; i < num_platforms; i++)
    {
      status = clGetPlatformInfo(list_platforms[i], CL_PLATFORM_VENDOR, sizeof(vendor_buffer), vendor_buffer, NULL);
      check_error(status, "Get platform info");

      if (strcmp(vendor_buffer, amd_vendor) == 0)
      {
        *platform = list_platforms[i];
        break;
      }
    }

    free(list_platforms);
  } 
}



int set_device_and_context(cl_platform_id *platform, cl_device_id *device, cl_context *context)
{
  int status;

  uint          num_devices;
  cl_device_id *list_devices;

  uint           device_number = 0; /* Get first device it sees */
  cl_device_type device_type   = CL_DEVICE_TYPE_GPU;

  status = clGetDeviceIDs(*platform, device_type, 0, NULL, &num_devices);
  check_error(status, "Get device ID");

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

  get_amd_platform(&platform);

  set_device_and_context(&platform, &device, context);

  *list_queues = malloc(num_queues * sizeof(cl_command_queue));

 
  for (i = 0; i < num_queues; i++)
  {
    (*list_queues)[i] = clCreateCommandQueue(*context, device, 0, &status);
    check_error(status, "Create command queue");
  }
}



/**
 * Attaches the received parameter as a kernel argument.-
 *
 */
void set_kernel_arg (cl_kernel *kernel, 
                     cl_uint arg_index, 
                     size_t arg_size, 
                     const void *arg_value)
{
    cl_int status;
    char msg [1024];

    snprintf (msg, 1024, "Set %d kernel argument", arg_index);

    status = clSetKernelArg (*kernel, 
                             arg_index,
                             arg_size,
                             arg_value);
    check_error (status, msg);
}



/**
 * Marks a kernel parameter as local memory.-
 *
 */
void set_local_mem (cl_kernel *kernel, cl_uint arg_index, size_t arg_size)
{
    cl_int status;
    char msg [1024];

    snprintf (msg, 
              1024, 
              "Set %d kernel argument as local memory, with size %ld",
              arg_index,
              arg_size);

    status = clSetKernelArg (*kernel, 
                             arg_index,
                             arg_size,
                             NULL);
    check_error (status, msg);
}



/**
 * Creates an OpenCL buffer, returning a pointer to it in the last parameter.-
 *
 */
void create_buffer (cl_context *context, cl_mem_flags flags, size_t size, void *host_ptr, cl_mem *dev_ptr)
{
    cl_int status;
    *dev_ptr = clCreateBuffer (*context,
                               flags,
                               size,
                               host_ptr,
                               &status);
    if (*dev_ptr == NULL)
        check_error (status, "Create buffer");
}



/**
 * Reads an OpenCL buffer from the device.-
 *
 */
void read_buffer (cl_command_queue *queue, 
                  cl_mem *dev_ptr, 
                  cl_bool blocking_read, 
                  size_t size, 
                  void *host_ptr)
{
    cl_int status;
    status = clEnqueueReadBuffer (*queue,
                                  *dev_ptr,
                                   blocking_read,
                                   0,
                                   size,
                                   host_ptr,
                                   0,
                                   NULL,
                                   NULL);
     check_error (status, "Read buffer");
}



/**
 * Writes an OpenCL buffer to the device.-
 *
 */
void write_buffer (cl_command_queue *queue, 
                   cl_mem *dev_ptr, 
                   cl_bool blocking_write, 
                   size_t size, 
                   const void *host_ptr)
{
    cl_int status;
    status = clEnqueueWriteBuffer (*queue,
                                   *dev_ptr,
                                   blocking_write,
                                   0,
                                   size,
                                   host_ptr,
                                   0,
                                   NULL,
                                   NULL);
     check_error (status, "Write buffer");
}



/**
 * Executes the specified kernel.-
 *
 */
void run_kernel (cl_command_queue *queue, cl_kernel *kernel)
{
    cl_int status;
    status = clEnqueueTask (*queue,
  	                        *kernel,
                            0,
                            NULL,
                            NULL);
     check_error (status, "Run kernel");
}



void release_opencl (int num_queues, 
                     cl_command_queue **list_queues, 
                     cl_context *context)
{
  int i;

  /* Flushes and finishes all queues in the list of queues  */
  for (i = 0; i < num_queues; i++)
  {
    clFinish((*list_queues)[i]);
  }

  clReleaseContext(*context);

  /* Free memory used by list of pointers */
  free(*list_queues);
}



/**
 * Builds a kernel for a given device and context from a source file.-
 *
 */
void build_kernel_from_file (cl_context *context, 
                             char *kernel_name, 
                             const char *file_name, 
                             cl_kernel *kernel)
{
    int bytes_read;
    char *source = (char *) calloc (1024000,
                                    sizeof (char));
    const char **kernel_source = (const char **) &(source[0]);

    bytes_read = read_kernel_file (file_name, 
                                   source);
    if (bytes_read < 1)
    {
        fprintf (stderr, 
                 "ERROR Could load kernel source from <%s>",
                 file_name);
        exit (1);
    }
    build_kernel (context,
                  kernel_name,
                  kernel_source,
                  kernel);
    free (source);
}



void build_kernel (cl_context *context, 
                   char *kernel_name, 
                   const char **kernel_source, 
                   cl_kernel *kernel)
{
  int          status;
  cl_device_id device;
  uint         num_devices = 1;

  get_device_from_context(context, &device);

  uint num_strings_in_kernel_source = 1;

  cl_program program = clCreateProgramWithSource(*context, num_strings_in_kernel_source, kernel_source, NULL, &status);
  check_error(status, "Create Program With Source");

  int build_status   = clBuildProgram(program, num_devices, &device, NULL, NULL, NULL);
  check_error(status, "Build program");

  if (build_status != CL_SUCCESS)
  {
    int error_buffer_size = 4000;
    char error_buffer[error_buffer_size];

    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, error_buffer_size, error_buffer, NULL);

    printf("OpenCL Build log:\n\n%s.\n", error_buffer);
  }

  *kernel = clCreateKernel(program, kernel_name, NULL);
}


/* EOF */
