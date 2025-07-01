#!/usr/bin/env python3
"""
Logger configuration module for Orv-Mol.
Provides a loguru-compatible interface using standard Python logging.
"""

import logging
import sys
from pathlib import Path
from datetime import datetime

class LoguruCompatibleLogger:
    """Logger class that provides loguru-like interface using standard logging."""
    
    def __init__(self):
        self.logger = logging.getLogger('orv_mol')
        self.logger.setLevel(logging.DEBUG)
        
        # Remove existing handlers to avoid duplicates
        self.logger.handlers.clear()
        
        # Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter(
            '%(asctime)s | %(levelname)8s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        console_handler.setFormatter(console_formatter)
        self.logger.addHandler(console_handler)
        
        # File handler (if logs directory exists)
        logs_dir = Path('logs')
        if logs_dir.exists() or self._create_logs_dir():
            log_file = logs_dir / f'orv_mol_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(logging.DEBUG)
            file_formatter = logging.Formatter(
                '%(asctime)s | %(levelname)8s | %(name)s:%(lineno)d | %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
            file_handler.setFormatter(file_formatter)
            self.logger.addHandler(file_handler)
    
    def _create_logs_dir(self):
        """Create logs directory if possible."""
        try:
            Path('logs').mkdir(exist_ok=True)
            return True
        except Exception:
            return False
    
    def info(self, message):
        """Log info message."""
        self.logger.info(message)
    
    def debug(self, message):
        """Log debug message."""
        self.logger.debug(message)
    
    def warning(self, message):
        """Log warning message."""
        self.logger.warning(message)
    
    def error(self, message):
        """Log error message."""
        self.logger.error(message)
    
    def success(self, message):
        """Log success message (mapped to info with special prefix)."""
        self.logger.info(f"SUCCESS: {message}")
    
    def critical(self, message):
        """Log critical message."""
        self.logger.critical(message)

def setup_logging(log_level: str = "INFO", log_file: str = None):
    """Setup logging configuration with specified level and optional file output."""
    # Convert string log level to logging constant
    level_map = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
        "CRITICAL": logging.CRITICAL
    }
    
    level = level_map.get(log_level.upper(), logging.INFO)
    
    # Get the logger instance
    logger_obj = logging.getLogger('orv_mol')
    logger_obj.setLevel(level)
    
    # Clear existing handlers
    logger_obj.handlers.clear()
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_formatter = logging.Formatter(
        '%(asctime)s | %(levelname)8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_handler.setFormatter(console_formatter)
    logger_obj.addHandler(console_handler)
    
    # File handler (custom or auto-generated)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter(
            '%(asctime)s | %(levelname)8s | %(name)s:%(lineno)d | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(file_formatter)
        logger_obj.addHandler(file_handler)
    else:
        # Auto file logging in logs directory
        logs_dir = Path('logs')
        try:
            logs_dir.mkdir(exist_ok=True)
            log_file_auto = logs_dir / f'orv_mol_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
            file_handler = logging.FileHandler(log_file_auto)
            file_handler.setLevel(logging.DEBUG)
            file_formatter = logging.Formatter(
                '%(asctime)s | %(levelname)8s | %(name)s:%(lineno)d | %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
            file_handler.setFormatter(file_formatter)
            logger_obj.addHandler(file_handler)
        except Exception:
            pass  # File logging is optional

# Create global logger instance
logger = LoguruCompatibleLogger() 