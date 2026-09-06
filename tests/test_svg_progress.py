import io
import unittest
from unittest import mock

from svg_progress import SvgProgressBar


class TerminalStream(io.StringIO):
    def isatty(self):
        return True


def visible_line(text):
    line, column = [], 0
    for char in text:
        if char == "\r":
            column = 0
        else:
            if column == len(line):
                line.append(char)
            else:
                line[column] = char
            column += 1
    return "".join(line).strip()


class SvgProgressTests(unittest.TestCase):
    def test_bar_is_erased_on_completion(self):
        stream = TerminalStream()
        with SvgProgressBar(stream) as progress:
            progress.update(1, 3)
            self.assertIn("1/3 submitted pages", visible_line(stream.getvalue()))
            progress.update(3, 3)
            self.assertEqual(visible_line(stream.getvalue()), "")
        self.assertNotIn("\n", stream.getvalue())

    def test_bar_is_erased_on_exception(self):
        stream = TerminalStream()
        with self.assertRaises(RuntimeError):
            with SvgProgressBar(stream) as progress:
                progress.update(0, 4)
                raise RuntimeError("failed")
        self.assertEqual(visible_line(stream.getvalue()), "")

    def test_redirected_output_has_no_progress_lines(self):
        stream = io.StringIO()
        with SvgProgressBar(stream) as progress:
            progress.update(0, 4)
            progress.update(4, 4)
        self.assertEqual(stream.getvalue(), "")

    def test_step_log_is_printed_after_erasing_bar(self):
        stream = TerminalStream()
        log = mock.Mock(
            side_effect=lambda result: self.assertEqual(
                visible_line(stream.getvalue()), "",
            )
        )
        with SvgProgressBar(stream, print_step=log) as progress:
            progress.update(1, 2)
            progress.log_step("step")
        log.assert_called_once_with("step")
