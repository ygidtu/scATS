#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created by Zhang
"""
from rich.progress import (BarColumn, Progress, TextColumn, TimeElapsedColumn,
                            TimeRemainingColumn, TransferSpeedColumn)


def custom_progress(io: bool = False):
    if io:
        return Progress(
            "[progress.description]{task.description}",
            BarColumn(),
            "[progress.percentage]{task.percentage:>3.0f}%",
            TextColumn("| Remaining:"),
            TimeRemainingColumn(),
            TextColumn("| Speed: "),
            TransferSpeedColumn()
        )
    return Progress(
        "[progress.description]{task.description}",
        BarColumn(),
        "[progress.percentage]{task.percentage:>3.0f}% ({task.completed}/{task.total})",
        TextColumn("| Elapsed:"),
        TimeElapsedColumn(),
        TextColumn("| Remaining:"),
        TimeRemainingColumn(),
    )


if __name__ == '__main__':
    pass
