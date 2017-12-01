

#define FE_MESSAGE_BUFFER_SIZE 1024
/**
 * Exception class with printf like text formatting
 * capabilities, using variable argument list.
 */
class FException {
private:
	char error_msg[FE_MESSAGE_BUFFER_SIZE];
	unsigned int error_code;

protected:

public:
	FException(char* error_msg, ...);
	FException(unsigned int error_code, char* error_msg, ...);

	void Print();
	unsigned int getErrorCode();
	char *getMessage();
};
