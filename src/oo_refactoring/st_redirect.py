# Acknowledgement:  Bela Schaum for original code
# Reference:        https://github.com/streamlit/streamlit/issues/268

import contextlib
import io

import streamlit as st


class StreamlitRedirect:
    """
    A class to redirect stdout or stderr to a Streamlit container.
    """

    class StreamBuffer(io.StringIO):
        def __init__(self, output_func, max_buffer, buffer_separator):
            super().__init__()
            self.output_func = output_func
            self.max_buffer = max_buffer
            self.buffer_separator = buffer_separator

        def write(self, s):
            if self.max_buffer:
                current_length = super().tell()
                if current_length + len(s) > self.max_buffer:
                    # Truncate buffer if it exceeds the max_buffer size
                    super().seek(0)
                    super().write(s)
                    super().truncate()
                else:
                    super().write(s)
            else:
                super().write(s)
            self.output_func(self.getvalue())
            return len(s)

    def __init__(self, to=None, to_file=None, max_buffer=None, buffer_separator="\n"):
        self.stream_buffer = self.StreamBuffer(
            self.write_to_streamlit, max_buffer, buffer_separator
        )
        self.to = to or st
        self.to_file = to_file
        self.container = None

    def __enter__(self):
        if self.container is not None:
            raise RuntimeError("This redirector is already in use")
        self.container = self.to.empty()
        self.context_manager = contextlib.redirect_stdout(self.stream_buffer)
        self.context_manager.__enter__()
        return self.stream_buffer

    def __exit__(self, exc_type, exc_value, traceback):
        self.context_manager.__exit__(exc_type, exc_value, traceback)
        self.context_manager = None
        if self.to_file:
            with open(self.to_file, "w") as file:
                file.write(self.stream_buffer.getvalue())
        self.container = None
        self.stream_buffer = self.stream_buffer.__class__(
            self.stream_buffer.output_func,
            self.stream_buffer.max_buffer,
            self.stream_buffer.buffer_separator,
        )

    def write_to_streamlit(self, data):
        self.container.code(data)  # Using 'code' for monospaced font


# Usage example:
# with StreamlitRedirect(to=st.sidebar, to_file="output.log"):
#     print("This will be displayed in the sidebar and written to 'output.log'")
