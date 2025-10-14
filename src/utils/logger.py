import os
import sys
import psutil
from loguru import logger

def setup_logger(app_name="phd_thesis"):
    """Configure the application logger.
    
    Args:
        app_name (str): Name to use for the log file (default: "phd_thesis")
    """
    logger.remove()  # Remove default handler
    
    # Create logs directory if it doesn't exist
    log_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "logs")
    os.makedirs(log_path, exist_ok=True)
    
    # Add file handler
    logger.add(
        os.path.join(log_path, f"{app_name}.log"),
        rotation="100 MB",
        retention="10 days",
        level="INFO",
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>"
    )
    
    # Add console handler
    logger.add(
        sys.stderr,
        level="INFO",
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>"
    )
    
    return logger

def get_memory_usage():
    """Get current memory usage of the process in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024  # Convert to MB

def log_memory_usage(context=""):
    """Log current memory usage with optional context."""
    mem_usage = get_memory_usage()
    logger.info(f"Memory usage{' ' + context if context else ''}: {mem_usage:.2f} MB")
    return mem_usage

def monitor_memory(func):
    """Decorator to monitor memory usage before and after a function call."""
    from functools import wraps
    
    @wraps(func)
    def wrapper(*args, **kwargs):
        func_name = func.__name__
        start_mem = log_memory_usage(f"before {func_name}")
        
        try:
            result = func(*args, **kwargs)
            end_mem = log_memory_usage(f"after {func_name}")
            delta_mem = end_mem - start_mem
            
            if delta_mem > 100:  # Log warning if function used more than 100MB
                logger.warning(f"High memory usage in {func_name}: {delta_mem:.2f} MB")
            else:
                logger.info(f"Memory delta for {func_name}: {delta_mem:.2f} MB")
                
            return result
        except Exception as e:
            logger.error(f"Error in {func_name}: {str(e)}")
            raise
    
    return wrapper
